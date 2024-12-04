#include "force_chain.hpp"
#include <image/nrrd_wrapper.hpp>
#include <fstream>
#include <math/stat.hpp>
#include <teem/hest.h>

sparse_matrix<float>    forces;
std::vector<vec3>       coords;
double                  max_angle;
const double            particle_radius = 0.05;

char* coord_f, *force_f, *stress_f, *chain_f;
int N;
float a;
void initialize(int argc, const char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "c",      "coordinates",          airTypeString,  1, 1, &coord_f,             NULL,       "input coordinates file name");
    hestOptAdd(&hopt, "f",      "forces",               airTypeString,  1, 1, &force_f,             NULL,       "input forces file name");
    hestOptAdd(&hopt, "o",      "chain output",         airTypeString,  1, 1, &chain_f,             NULL,       "output force chain file name");
    hestOptAdd(&hopt, "s",      "stress output",        airTypeString,  1, 1, &stress_f,            NULL,       "output stress file");
    hestOptAdd(&hopt, "n",      "# valid points",       airTypeInt,     1, 1, &N,                   NULL,       "number of valid particles");
    hestOptAdd(&hopt, "a",      "max angle",            airTypeFloat,   0, 1, &a,                   "60",       "max angular tolerance");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Extract force chains in granular particle assembly and export principal stress tensors",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

void get_neighbors(std::vector<unsigned int>& neighbors, unsigned int i)
{
    std::vector<std::pair<unsigned int, float> > row;
    forces.row_entries(row, i);
    neighbors.resize(row.size());
    for (int l = 0 ; l < row.size() ; ++l) {
        neighbors[l] = row[l].first;
    }
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& values)
{
    std::pair<double, double> mv = spurt::meanvariance(values);
    T min = *std::min_element(values.begin(), values.end());
    T max = *std::max_element(values.begin(), values.end());
    os << " \tmin=" << min << ", \tmax=" << max << ", \tmean=" << mv.first
       << ", \tstd dev=" << sqrt(mv.second);
   return os;
}

std::pair<double, double> contact_pressures(const nvis::fvec3 contact,
        const nvis::fvec3& center, float force)
{
    nvis::fvec3 radial_dir(contact[0], contact[1], 0);
    radial_dir /= nvis::norm(radial_dir);
    nvis::fvec3 vec_f = contact - center;
    vec_f *= force / nvis::norm(vec_f);
    
    double radial_p = fabs(nvis::inner(vec_f, radial_dir));
    double axial_p = fabs(vec_f[2]);
    
    return std::make_pair(radial_p, axial_p);
}

std::pair<double, double> particle_pressures(unsigned int idx)
{
    double axial = 0;
    double radial = 0;
    
    std::vector<std::pair<unsigned int, float> > row;
    forces.row_entries(row, idx);
    for (unsigned int i = 0 ; i < row.size() ; ++i) {
        nvis::fvec3 contact = 0.5 * (coords[idx], coords[row[i].first]);
        std::pair<double, double> p = contact_pressures(contact, coords[idx],
                                      row[i].second);
        radial += p.first;
        axial += p.second;
    }
    
    return std::make_pair(radial, axial);
}

int main(int argc, const char* argv[])
{
    initialize(argc, argv);
    
    chain_type
    current_chain,          // chain being currently computed
    current_front;          // set of positions from which to grow the chain in next iteration
    
    std::vector<chain_type> all_chains;
    
    // NB:  this value corresponds to the average stress value
    double sigma_mean = 0;
    // NB: these values corresponds to radial and axial pressure
    double total_radial_p = 0;
    double total_axial_p = 0;
    double chain_tangent_radial_p = 0;
    double chain_tangent_axial_p = 0;
    double chain_contact_radial_p = 0;
    double chain_contact_axial_p = 0;
    std::vector<particle> seed_particles; // all admissible particles from which to start a chain
    
    Nrrd* contacts, *pos;
    pos = spurt::readNrrd(coord_f);
    contacts = spurt::readNrrd(force_f);
    max_angle = M_PI * a / 180.;
    
    // convert nrrds to internal data structure
    coords.resize(N);
    std::cerr << "coordinate file contains " << pos->axis[1].size << " particle coordinates\n";
    float* data = (float*)pos->data;
    for (int i = 0 ; i < N ; ++i) {
        coords[i][0] = data[3*i  ];
        coords[i][1] = data[3*i+1];
        coords[i][2] = data[3*i+2];
    }
    
    std::vector < std::list<std::pair<unsigned int, float> > > all_forces(N);
    data = (float*)contacts->data;
    int n_contacts = contacts->axis[1].size;
    for (int l = 0 ; l < n_contacts ; ++l) {
        unsigned int i  = data[3*l  ] - 1;
        unsigned int j  = data[3*l+1] - 1;
        float f         = data[3*l+2];
        
        nvis::fvec3 contact;
        if (j < N) {
            all_forces[i].push_back(std::pair<unsigned int, float>(j, f));
            all_forces[j].push_back(std::pair<unsigned int, float>(i, f));
            
            nvis::fvec3 contact = 0.5 * (coords[i] + coords[j]);
        } else {
            j -= N;
            switch (j) {
                case 0: {
                    // std::cerr << "bottom face\n";
                    // bottom face
                    contact = coords[i];
                    contact[2] -= 2.*particle_radius;
                    coords.push_back(contact);
                    unsigned int k = coords.size() - 1;
                    all_forces[i].push_back(std::pair<unsigned int, float>(k, f));
                    break;
                }
                case 1: {
                    // std::cerr << "side face\n";
                    // side face
                    contact = coords[i];
                    nvis::fvec3 radial(coords[i][0], coords[i][1], 0);
                    radial *= particle_radius / nvis::norm(radial);
                    contact += 2.*radial;
                    coords.push_back(contact);
                    unsigned int k = coords.size() - 1;
                    all_forces[i].push_back(std::pair<unsigned int, float>(k, f));
                    break;
                }
                case 2: {
                    // std::cerr << "top face\n";
                    // top face
                    contact = coords[i];
                    contact[2] += 2.*particle_radius;
                    coords.push_back(contact);
                    unsigned int k = coords.size() - 1;
                    all_forces[i].push_back(std::pair<unsigned int, float>(k, f));
                    break;
                }
                default: {
                    std::cerr << "unknown boundary case\n";
                    break;
                }
            }
            
            std::pair<double, double> p = contact_pressures(contact, coords[i], f);
            total_radial_p += p.first;
            total_axial_p += p.second;
        }
    }
    std::cerr << coords.size() - N << " boundary particles have been added\n";
    all_forces.resize(coords.size());
    
    std::cerr << "total axial pressure (not normalized) = " << total_axial_p << '\n';
    std::cerr << "total radial pressure (not normalized) = " << total_radial_p << '\n';
    
    forces.initialize(all_forces);
    std::cerr << "force matrix initialized\n";
    
    // compute mean stress value
    sigma_mean = 0;
    std::vector<particle> all_particles(N);
    for (int i = 0 ; i < N ; ++i) {
        // if (!(i % 1)) {
        //  std::cerr << "computing stress tensor for particle " << i << " / " << N
        //            << std::endl;
        // }
        sym_mat3 stress = compute_stress(i);
        nvis::vec3 dir;
        double val, cl;
        eigenanalysis(dir, val, cl, stress);
        
        all_particles[i].index = i;
        all_particles[i].stress_direction = dir;
        all_particles[i].principal_stress = val;
        all_particles[i].stress = stress;
        all_particles[i].cl = cl;
        
        sigma_mean += val;
    }
    sigma_mean /= (double)N;
    
    // export stress value
    {
        float* stress_data = (float*)calloc(6 * N, sizeof(float));
        for (int i = 0 ; i < N ; ++i) {
            for (int j = 0 ; j < 6 ; ++j) {
                stress_data[6*i+j] = all_particles[i].stress[j];
            }
        }
        std::vector<size_t> dims(2);
        dims[0] = 6;
        dims[1] = N;
        std::vector<double> spc(2);
        spc[0] = spc[1] = airNaN();
        spurt::writeNrrdFromContainers(stress_data, stress_f, /*nrrdTypeFloat,*/ dims, spc);
        std::cerr << "stress data exported\n";
    }
    
    seed_particles.clear();
    for (int i = 0 ; i < N ; ++i) {
        if (all_particles[i].principal_stress >= sigma_mean) {
            seed_particles.push_back(all_particles[i]);
        }
    }
    
    std::cout << "mean principal stress = " << sigma_mean << ", "
              << seed_particles.size() << " particles over threshold ("
              << 100*seed_particles.size() / all_particles.size() << "%)\n";
              
    std::sort(seed_particles.begin(), seed_particles.end(), Lt_stress());
    
    std::set<unsigned int> marked;
    std::vector<float> chain_total_stress;
    std::vector<float> chain_avg_stress;
    std::vector<float> chain_length;
    std::vector<float> chain_total_anisotropy;
    std::vector<float> chain_avg_anisotropy;
    
    // loop over all identified seeds
    for (int i = 0 ; i < seed_particles.size() ; ++i) {
    
        // std::cout << "processing seed " << i << " from " << seed_particles.size() << std::endl;
        
        const particle& p = seed_particles[i];
        if (marked.find(p.index) != marked.end()) {
            continue;
        }
        marked.insert(p.index);
        
        std::pair<double, double> pr = particle_pressures(p.index);
        chain_contact_radial_p = pr.first;
        chain_contact_axial_p = pr.second;
        
        int s = p.index;
        
        current_chain.clear();
        current_chain.insert(p);
        current_front.clear();
        
        std::vector<unsigned int> neighbors;
        get_neighbors(neighbors, s);
        vec3 sigma = p.stress_direction;
        
        // identify neighbor(s) belonging to force chain. they will form the initial front
        for (int i = 0 ; i < neighbors.size() ; ++i) {
            unsigned int n = neighbors[i];
            if (marked.find(n) == marked.end()) {
                // particle is not part of another force chain
                if (is_aligned(sigma, s, n, true) &&
                        (n >= N || is_aligned(all_particles[n].stress_direction, n, s, true))) {
                    // these particles are aligned with their respective stress direction
                    if (n < N) {
                        all_particles[n].connect_to = s;
                        marked.insert(n);
                        
                        if (all_particles[n].principal_stress >= sigma_mean) {
                            current_front.insert(all_particles[n]);
                        }
                        
                        current_chain.insert(all_particles[n]);
                    }
                }
            }
        }
        
        // extend chain from particles in front
        while (!current_front.empty()) {
            particle pic = *current_front.begin();
            current_front.erase(current_front.begin()); // remove current particle from front
            int id = pic.index;
            vec3 incoming_dir = radius(pic.connect_to, id); // direction from which we are coming
            sigma = pic.stress_direction;
            
            // compute radial and axial pressure
            nvis::fvec3 contact = 0.5 * (coords[pic.connect_to] + coords[pic.index]);
            std::pair<double, double> pr =
                contact_pressures(contact, coords[pic.index],
                                  forces(pic.index, pic.connect_to));
            chain_contact_radial_p += pr.first;
            chain_contact_axial_p += pr.second;
            
            nvis::fvec3 radial_dir(contact[0], contact[1], 0);
            radial_dir /= nvis::norm(radial_dir);
            
            // ensure that sigma is pointing in correct direction. if not, fix it.
            sigma *= (nvis::inner(sigma, incoming_dir) < 0 ? -1 : 1);
            
            // update the front with neighbors of current particle
            get_neighbors(neighbors, id);
            for (int i = 0 ; i < neighbors.size() ; ++i) {
                int n = neighbors[i];
                if (marked.find(n) == marked.end()) {
                    // particle is not part of another force chain
                    if (is_aligned(sigma, id, n, false) &&
                            // we are only interested in forward direction here
                            (n >= N || is_aligned(all_particles[n].stress_direction, n, id, true))) {
                        if (n < N) {
                            all_particles[n].connect_to = id;
                            marked.insert(n);
                            
                            if (all_particles[n].principal_stress >= sigma_mean) {
                                current_front.insert(all_particles[n]);
                            }
                            
                            current_chain.insert(all_particles[n]);
                        }
                    }
                }
            }
        }
        
        // did we find a true force chain?
        if (current_chain.size() >= 3) {
            all_chains.push_back(current_chain);
            
            // some basic statistics
            chain_total_stress.push_back(total_stress(current_chain));
            chain_avg_stress.push_back(average_stress(current_chain));
            chain_length.push_back(current_chain.size());
            chain_total_anisotropy.push_back(total_anisotropy(current_chain));
            chain_avg_anisotropy.push_back(average_anisotropy(current_chain));
        }
        // NB: pairs and singletons are tagged but not included as chains
    }
    
    std::cerr << "chains statistics\n";
    std::cerr << "total stress: " << chain_total_stress << '\n'
              << "average stress: " << chain_avg_stress << '\n'
              << "length: " << chain_length << '\n'
              << "total anisotropy: " << chain_total_anisotropy << '\n'
              << "average anisotropy: " << chain_avg_anisotropy << '\n';
              
    std::cerr << "\n\n"
              << "total radial pressure: " << total_radial_p << '\n'
              << "total axial pressure: " << total_axial_p << '\n'
              << "total pressure ratio: " << 100.*total_radial_p / total_axial_p << "%\n"
              << "chain contact radial pressure: " << chain_contact_radial_p << '\n'
              << "chain contact axial pressure: " << chain_contact_axial_p << '\n'
              << "contact chain pressure ratio: " << 100.*chain_contact_radial_p / chain_contact_axial_p << "%\n"
              << "chain tangent radial pressure: " << chain_tangent_radial_p << '\n'
              << "chain tangent axial pressure: " << chain_tangent_axial_p << '\n'
              << "tangent chain pressure ratio: " << 100.*chain_tangent_radial_p / chain_tangent_axial_p << "%\n";
              
              
    std::cout << "saving " << all_chains.size() << " to file\n";
    std::fstream output(chain_f, std::ios::out);
    output << "# Format description:\n"
           << "# number_of_particles_in_chain (unsigned int)\n"
           << "# particle_index (unsigned int) coordinates (3 floats) stress_value (float)...\n"
           << "# ...stress_direction (3 floats) stress_anisotropy (float)...\n"
           << "# ...previous_point_in_chain (unsigned int)\n";
    for (unsigned int i = 0 ; i < all_chains.size() ; ++i) {
        std::set<particle>::const_iterator iter;
        output << all_chains[i].size() << '\n';
        for (iter = all_chains[i].begin() ; iter != all_chains[i].end() ; ++iter) {
            const particle& p = *iter;
            output << p.index << " "
                   << coords[p.index][0] << " "
                   << coords[p.index][1] << " "
                   << coords[p.index][2] << " "
                   << p.principal_stress << " "
                   << p.stress_direction[0] << " "
                   << p.stress_direction[1] << " "
                   << p.stress_direction[2] << " "
                   << p.cl << " "
                   << p.connect_to << '\n';
        }
    }
    output.close();
}
















































































