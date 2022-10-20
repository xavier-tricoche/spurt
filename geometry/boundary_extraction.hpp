#ifndef __XAVIER_GEOMETRY_BOUNDARY_EXTRACTION_HPP__
#define __XAVIER_GEOMETRY_BOUNDARY_EXTRACTION_HPP__

#include <format/DLRreader.hpp>
#include <list>
#include <map>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <misc/progress.hpp>
#include <math/fixed_vector.hpp>
#include <image/nrrd_wrapper.hpp>

namespace spurt {

template<typename T>
struct subvector {
    typedef typename std::vector<T> container_type;
    typedef typename container_type::value_type value_type;
    typedef typename container_type::const_iterator const_iterator;

    subvector(const_iterator begin, const_iterator end);
    subvector(const container_type& c, size_t first, size_t length);

    const value_type& operator[](int i) const ;

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator m_begin, m_end;
};

class mesh {
public:
    typedef DLRreader::cell_type cell_kind;
    typedef long int index_type;
    static constexpr index_type invalid_index = std::numeric_limits<index_type>::max();
    typedef nvis::fixed_vector<float, 3> vertex_type;
    typedef nvis::fixed_vector<float, 3> vec3_type;
    typedef subvector<index_type> subvector_type;
    typedef std::pair<cell_kind, const subvector_type> cell_type;

    class face : public std::array<index_type, 4> {
    private:
        bool is_valid() const;
        bool is_triangle(bool check=false) const;
        bool is_quad(bool check=false) const;

    public:
        typedef std::array<index_type, 4> base_type;
        typedef base_type::iterator iterator;
        typedef base_type::const_iterator const_iterator;

        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;

        size_t size() const;

        face();
        face(const face& other);
        face(index_type i, index_type j, index_type k, int id, bool sort=true);
        face(index_type i, index_type j, index_type k, index_type l, int id, bool sort=true) ;
        void sort();
        bool operator<(const face& other) const;
        int face_id() const;
        int& face_id();

        int m_id;
    };

    mesh() {}
    mesh(const std::vector<vertex_type>& vertices, const std::vector<long>& cell_indices, const std::vector<std::pair<DLRreader::cell_type, long>>& cell_offsets)
        : m_vertices(vertices), m_cell_indices(cell_indices), m_cell_offsets(cell_offsets), m_ncells(cell_offsets.size()-1) {}

    void load_mesh(const std::string& mesh_name);
    void load_data(const std::string& data_name);

    cell_type get_cell(const index_type& id) const;

    static const std::vector<face>& reference_faces(cell_kind kind);

    void get_faces(const cell_type& cell, std::vector<face>& faces, bool sort=true);

    void get_actual_face(const cell_type& cell, int id, face& f);
    void extract_boundary(bool verbose=false);
    std::vector<face> get_boundary();
    std::vector<face> get_offset();
    std::vector<vec3_type> get_shear_stress(bool verbose=false);
    const std::vector<vertex_type>& get_vertices() const;
    const std::vector<vec3_type>& get_velocity() const;
    const std::vector<vec3_type>& get_normals() const;
    void save_boundary(const std::string& basename, bool verbose=false);
    void load_boundary(const std::string& basename, bool verbose=false);

private:
    static const std::vector<face> _tetrahedron;
    static const std::vector<face> _pyramid;
    static const std::vector<face> _prism;
    static const std::vector<face> _hexahedron;
    std::vector<vertex_type> m_vertices; // 3D mesh vertices
    std::vector<vec3_type>   m_velocity; // velocity vector field
    std::vector<long>  m_cell_indices;   // 1D list of all cells's vertices
    std::vector<std::pair<DLRreader::cell_type, long>> m_cell_offsets; // per-cell access in cell_indices
    size_t m_ncells; // number of cells
    std::map<face, std::list<index_type>> m_face_to_cells; // face-to-cell incidence information
    std::vector<std::pair<face, index_type>> m_boundary_to_cell; // boundary to cells relationship
    std::vector<std::pair<face, index_type>> m_offset_to_boundary; // offset surface to boundary correspondence
    std::vector<vec3_type> m_normals; // computed boundary normals
};

namespace {

    struct index_wrapper {
        index_wrapper(mesh::index_type id) : m_id(id) {}
        mesh::index_type m_id;
    };
}

static std::ostream& operator<<(std::ostream& os, const index_wrapper id) {
    if (id.m_id == mesh::invalid_index) {
        os << "<invalid>";
    }
    else {
        os << id.m_id;
    }
    return os;
}

static std::ostream& operator<<(std::ostream& os, const mesh::face& f) {
    os << "[";
    if (f.size() == 3) os << "triangle: ";
    else os << "quad: ";
    os << f.face_id() << "; ";
    std::for_each(f.begin(), f.end(), [&](mesh::index_type i) {
        os << index_wrapper(i) << ", ";
    });
    os << "]";
    return os;
}

static std::ostream& operator<<(std::ostream& os, const mesh::subvector_type& sv) {
    os << "[";
    std::for_each(sv.begin(), sv.end(), [&](mesh::index_type i) {
        os << index_wrapper(i) << ", ";
    });
    os << "]";
    return os;
}

static std::ostream& operator<<(std::ostream& os, const mesh::cell_type c) {
    std::string cellname;
    switch (c.first) {
        case DLRreader::TETRAHEDRON: cellname="tetrahedron"; break;
        case DLRreader::HEXAHEDRON: cellname="hexahedron"; break;
        case DLRreader::PYRAMID: cellname="pyramid"; break;
        case DLRreader::PRISM: cellname="prism"; break;
        default: cellname="unknown cell";
    }
    os << "(" << cellname << "," << c.second << ")";
    return os;
}

const std::vector<mesh::face> mesh::_tetrahedron
    { { 0, 1, 2,      0,   false },
      { 1, 3, 2,      1,   false },
      { 0, 2, 3,      2,   false },
      { 0, 3, 1,      3,   false } };

const std::vector<mesh::face> mesh::_pyramid
    { { 0, 1, 2, 3,   0,   false },
      { 0, 4, 1,      1,   false },
      { 1, 4, 2,      2,   false },
      { 2, 4, 3,      3,   false },
      { 0, 4, 3,      4,   false } };

const std::vector<mesh::face> mesh::_prism
    { { 0, 1, 2,      0,   false },
      { 3, 4, 5,      1,   false },
      { 0, 3, 5, 2,   2,   false },
      { 1, 2, 5, 4,   3,   false },
      { 0, 3, 4, 1,   4,   false } };

const std::vector<mesh::face> mesh::_hexahedron
    { { 0, 1, 2, 3,   0,   false },
      { 4, 5, 6, 7,   1,   false },
      { 0, 1, 5, 4,   2,   false },
      { 1, 2, 6, 5,   3,   false },
      { 3, 2, 6, 7,   4,   false },
      { 0, 3, 7, 4,   5,   false } };

// --- subvector ---

template<typename T>
subvector<T>::subvector(typename subvector<T>::const_iterator begin, typename subvector<T>::const_iterator end)
    : m_begin(begin), m_end(end) {}

template<typename T>
subvector<T>::subvector(const typename subvector<T>::container_type& c, size_t first, size_t length)
    : m_begin(c.begin() + first), m_end(c.end() + first + length) {}

template<typename T>
const typename subvector<T>::value_type& subvector<T>::operator[](int i) const {
    return *(m_begin + i);
}

template<typename T>
typename subvector<T>::const_iterator subvector<T>::begin() const {
    return m_begin;
}

template<typename T>
typename subvector<T>::const_iterator subvector<T>::end() const {
    return m_end;
}

// --- face ---

bool mesh::face::is_valid() const {
    return (*this)[0]!=invalid_index &&
           (*this)[1]!=invalid_index &&
           (*this)[2]!=invalid_index;
}
bool mesh::face::is_triangle(bool check) const {
    return (!check || is_valid()) && (*this)[3]==invalid_index;
}
bool mesh::face::is_quad(bool check) const {
    return (!check || is_valid()) && (*this)[3]!=invalid_index;
}
mesh::face::iterator mesh::face::begin() { return base_type::begin(); }
mesh::face::const_iterator mesh::face::begin() const { return base_type::begin(); }
mesh::face::iterator mesh::face::end() { return is_quad() ? base_type::end() : begin()+3; }
mesh::face::const_iterator mesh::face::end() const { return is_quad() ? base_type::end() : begin()+3; }

size_t mesh::face::size() const { return is_quad() ? 4 : 3; }

mesh::face::face() : base_type({invalid_index, invalid_index, invalid_index, invalid_index}), m_id(-1) {}
mesh::face::face(const face& other) : base_type(static_cast<base_type>(other)), m_id(other.m_id) {}
mesh::face::face(index_type i, index_type j, index_type k, int id, bool sort)
    : base_type({i, j, k, invalid_index}), m_id(id) {
    if (sort) {
        std::sort(this->begin(), this->begin()+3);
    }
}
mesh::face::face(index_type i, index_type j, index_type k, index_type l, int id, bool sort)
    : base_type({i, j, k, l}), m_id(id) {
    if (sort) {
        std::sort(this->begin(), this->end());
    }
}
void mesh::face::sort() { std::sort(begin(), end()); }
bool mesh::face::operator<(const mesh::face& other) const {
    for (int i=0; i<size(); ++i) {
        if ((*this)[i] < other[i]) return true;
        else if (other[i] < (*this)[i]) return false;
    }
    return false; // both faces are equal
}

int mesh::face::face_id() const { return m_id; }
int& mesh::face::face_id() { return m_id; }

// --- mesh ----

void mesh::load_mesh(const std::string& mesh_name) {
    spurt::DLRreader reader(mesh_name, "");
    reader.read_mesh(false, m_vertices, m_cell_indices, m_cell_offsets);
    m_ncells = m_cell_offsets.size()-1; // last entry is not an actual cell
}

void mesh::load_data(const std::string& data_name) {
    spurt::DLRreader reader("", data_name);
    reader.read_vector_data("velocity", m_velocity);
}

const std::vector<mesh::vertex_type>& mesh::get_vertices() const {
    if (m_vertices.empty()) {
        throw std::runtime_error("Error: No vertices available");
    }
    return m_vertices;
}

const std::vector<mesh::vec3_type>& mesh::get_velocity() const {
    if (m_velocity.empty()) {
        throw std::runtime_error("Error: No velocity data available. Did you forget to load the data?");
    }
    return m_velocity;
}

mesh::cell_type mesh::get_cell(const index_type& id) const {
    DLRreader::cell_type kind = m_cell_offsets[id].first;
    auto begin = m_cell_indices.begin() + m_cell_offsets[id].second;
    auto end = m_cell_indices.begin() + m_cell_offsets[id+1].second;
    return std::make_pair(kind, subvector_type(begin, end));
}

const std::vector<mesh::face>& mesh::reference_faces(mesh::cell_kind kind) {
    switch (kind) {
        case DLRreader::TETRAHEDRON: return _tetrahedron;
        case DLRreader::PYRAMID: return _pyramid;
        case DLRreader::PRISM: return _prism;
        case DLRreader::HEXAHEDRON: return _hexahedron;
        default: throw std::runtime_error("Unknown cell type" + std::to_string(kind));
    }
}

void mesh::get_faces(const mesh::cell_type& cell, std::vector<mesh::face>& faces, bool sort) {
    const std::vector<face>& ref = reference_faces(cell.first);
    faces.resize(ref.size());
    for (int n=0; n<ref.size(); ++n) {
        const face& rf = ref[n];
        for (int i=0; i<rf.size(); ++i) {
            faces[n][i] = cell.second[rf[i]];
        }
        faces[n].face_id() = rf.face_id();
        if (sort) faces[n].sort();
    }
}

void mesh::get_actual_face(const mesh::cell_type& cell, int id, mesh::face& f) {
    const std::vector<face>& ref = reference_faces(cell.first);
    const face& rf = ref[id];
    for (int i=0; i<rf.size(); ++i) {
        f[i] = cell.second[rf[i]];
    }
    f.face_id() = rf.face_id();
}

namespace {
    double norm(const mesh::vec3_type& n) {
        return nvis::norm(n);
    }

    mesh::vec3_type& normalize(mesh::vec3_type& v) {
        v /= nvis::norm(v);
        return v;
    }

    mesh::vec3_type& match(mesh::vec3_type& v, const mesh::vec3_type& r) {
        double d = nvis::inner(v, r);
        if (d < 0) {
            v *= -1;
        }
        return v;
    }

    bool is_zero(const mesh::vec3_type& v) {
        return v[0]==0 && v[1]==0 && v[2]==0;
    }

    mesh::vec3_type normal(const mesh::vertex_type& v0,
                           const mesh::vertex_type& v1,
                           const mesh::vertex_type& v2) {
        mesh::vec3_type v10 = v1-v0;
        mesh::vec3_type v20 = v2-v0;
        mesh::vec3_type n = nvis::cross(v10, v20);
        return normalize(n);
    }

    double dot(const mesh::vec3_type& v0, const mesh::vec3_type& v1) {
        return nvis::inner(v0, v1);
    }
}

void mesh::extract_boundary(bool verbose) {
    if (verbose) std::cout << "Entering mesh::extract_boundary()\n";
    m_boundary_to_cell.clear();
    m_offset_to_boundary.clear();
    m_face_to_cells.clear();
    spurt::ProgressDisplay progress;
    progress.fraction_on();
    progress.start(m_ncells, "Extracting all cell faces", 1000);
    for (index_type i=0; i<m_ncells; ++i) {
        progress.update(i+1);
        cell_type cell = get_cell(i);
        std::vector<face> faces;
        get_faces(cell, faces);
        std::for_each(faces.begin(), faces.end(), [&](const face& f) {
            auto iter = m_face_to_cells.find(f);
            if (iter == m_face_to_cells.end()) {
                m_face_to_cells[f] = std::list<long>();
            }
            m_face_to_cells[f].push_back(i);
        });
    } // loop over all cells
    progress.end();

    if (verbose)
        std::cout << "after processing all cells, there are "
              << m_face_to_cells.size() << " faces in this mesh\n";

    size_t ntri_prisms = 0;

    size_t faceid=0;
    progress.start(m_face_to_cells.size(), "Checking all faces", 1000);
    std::for_each(m_face_to_cells.begin(), m_face_to_cells.end(),
        [&](const std::pair<face, std::list<long>>& face_to_cell)
        {
            ++faceid;
            progress.update(faceid);
            if (face_to_cell.second.size() == 1) {
                // std::cout << "found a boundary face\n";
                const face& _f = face_to_cell.first;
                const index_type& cell_id = face_to_cell.second.front();
                // std::cout << "corresponding cell id is " << cell_id << '\n';
                cell_type cell = get_cell(cell_id);
                // if (cell.first == Prism && _f.face_id()>=0 && _f.face_id()<2) ++ntri_prisms;
                // std::cout << "cell is " << cell << '\n';
                // std::cout << "boundary face is " << std::flush;
                // std::cout << _f << '\n';
                face f;
                get_actual_face(cell, _f.face_id(), f);
                // std::cout << "actual face is " << f << '\n';
                m_boundary_to_cell.push_back(std::make_pair(f, cell_id));
                if (cell.first == DLRreader::PRISM && _f.face_id()>=0 && _f.face_id()<2) {
                    face of;
                    ++ntri_prisms;
                    get_actual_face(cell, _f.face_id() == 0 ? 1 : 0, of);
                    m_offset_to_boundary.push_back(std::make_pair(of, m_boundary_to_cell.size()-1));
                }
                else if (cell.first == DLRreader::HEXAHEDRON) {
                    int id=-1;
                    switch (_f.face_id()) {
                        case 0: id=1; break;
                        case 1: id=0; break;
                        case 2: id=4; break;
                        case 3: id=5; break;
                        case 4: id=2; break;
                        case 5: id=3; break;
                    }
                    if (id>=0) {
                        face of;
                        get_actual_face(cell, id, of);
                        m_offset_to_boundary.push_back(std::make_pair(of, m_boundary_to_cell.size()-1));
                    }
                }
            }
        }
    );
    progress.end();
    if (verbose) {
        std::cout << "There are " << ntri_prisms << " triangular faces of prisms on the boundary\n";
        std::cout << "There are " << m_offset_to_boundary.size() << " cells on 1-offset surface\n";
    }
    std::map<index_type, std::list<index_type>> vertex_to_cells;
    typedef std::map<index_type, std::list<index_type>>::value_type id_to_ids_type;

    std::vector<vec3_type> face_normals(m_boundary_to_cell.size());
    // compute connectivity information and cell-wise normals on the boundary
    progress.start(m_boundary_to_cell.size(), "Computing triangle-based normals", 1000);
    for (size_t i=0; i<m_boundary_to_cell.size(); ++i) {
        const face& f = m_boundary_to_cell[i].first;
        std::vector<vertex_type> verts(f.size());
        progress.update(i+1);
        for (int j=0; j<f.size(); ++j) {
            const index_type& ptid = f[j];
            auto iter = vertex_to_cells.find(ptid);
            if (iter == vertex_to_cells.end()) {
                vertex_to_cells[ptid] = std::list<index_type>();
            }
            vertex_to_cells[ptid].push_back(i);
            verts[j] = m_vertices[ptid];
        }

        if (f.size() == 3) {
            face_normals[i] = normal(verts[0], verts[1], verts[2]);
        }
        else {
            vec3_type n012 = normal(verts[0], verts[1], verts[2]);
            vec3_type n023 = normal(verts[0], verts[2], verts[3]);
            face_normals[i] = n012 + n023;
            normalize(face_normals[i]);
        }
    }
    progress.end();

    // compute boundary normals
    std::map<index_type, vec3_type> boundary_normals;
    typedef std::map<index_type, vec3_type>::value_type id2vec_type;
    size_t dcount=0;
    progress.start(vertex_to_cells.size(), "Compute vertex-based normals", 1000);
    std::for_each(vertex_to_cells.begin(), vertex_to_cells.end(),
        [&](const id_to_ids_type& i2is) {
            ++dcount;
            progress.update(dcount);
            const index_type& vid = i2is.first;
            const std::list<index_type>& cids = i2is.second;
            if (cids.empty()) {
                if (verbose) std::cerr << "Error: no face incident to boundary vertex\n";
                boundary_normals[vid] = vec3_type(0,0,0);
            }
            else {
                auto iter = cids.begin();
                vec3_type n0 = face_normals[*iter];
                vec3_type sum = n0;
                for (++iter; iter!=cids.end(); ++iter) {
                    vec3_type n1 = face_normals[*iter];
                    n1 = match(n1, n0);
                    sum += n1;
                }
                normalize(sum);
                boundary_normals[vid] = sum;
            }
        }
    );
    progress.end();

    m_normals.resize(m_vertices.size());
    std::fill(m_normals.begin(), m_normals.end(), vec3_type(0,0,0));
    std::for_each(boundary_normals.begin(), boundary_normals.end(), [&](const id2vec_type& id2v) {
       m_normals[id2v.first] = id2v.second;
    });


    if (verbose) std::cout << "Leaving mesh::extract_boundary()\n";

} // extract_boundary

std::vector<mesh::face> mesh::get_boundary() {
    if (m_boundary_to_cell.empty()) extract_boundary();
    std::vector<face> boundary(m_boundary_to_cell.size());
    for (size_t i=0; i<boundary.size(); ++i) {
        boundary[i]= m_boundary_to_cell[i].first;
    }
    return boundary;
}

std::vector<mesh::face> mesh::get_offset() {
    if (m_offset_to_boundary.empty()) extract_boundary();
    std::vector<face> offset(m_offset_to_boundary.size());
    for (size_t i=0; i<offset.size(); ++i) {
        offset[i] = m_offset_to_boundary[i].first;
    }
    return offset;
}

const std::vector<mesh::vec3_type>& mesh::get_normals() const {
    if (m_normals.empty()) throw std::runtime_error("No normals available");
    return m_normals;
}

std::vector<mesh::vec3_type> mesh::get_shear_stress(bool verbose) {
    if (verbose) std::cout << "\nentering mesh::get_shear_stress()\n";
    if (m_velocity.empty()) {
        throw std::runtime_error("Error: Unable to compute shear stress: no velocity data available");
    }
    if (m_boundary_to_cell.empty() || m_offset_to_boundary.empty() || m_normals.empty()) extract_boundary();

    spurt::ProgressDisplay progress;
    progress.fraction_on();
    std::map<index_type, vec3_type> boundary_velocity;
    std::vector<vec3_type> shear_stress(m_velocity.size());
    std::fill(shear_stress.begin(), shear_stress.end(), vec3_type(0,0,0));
    size_t nn=0;
    progress.start(m_offset_to_boundary.size(), "Computing shear stress vector field", 100);
    std::for_each(m_offset_to_boundary.begin(), m_offset_to_boundary.end(),
        [&](const std::pair<face, index_type>& f2f) {
            // std::cout << "m_offset_to_boundary: " << nn++ << ": " << f2f.first << " -> " << f2f.second << '\n';

            progress.update(++nn);

            const face& of = f2f.first;
            const face& bf = m_boundary_to_cell[f2f.second].first;
            for (int i=0; i<bf.size(); ++i) {
                auto iter = boundary_velocity.find(bf[i]);
                if (iter == boundary_velocity.end()) {
                    vec3_type v = m_velocity[of[i]]; // velocity at matching offset vertex
                    v -= dot(v, m_normals[bf[i]])*m_normals[bf[i]];
                    boundary_velocity[bf[i]] = v;
                    shear_stress[bf[i]] = v;
                }
            }
        }
    );
    progress.end();

    if (verbose) std::cout << "\nleaving mesh::get_shear_stress()\n";

    return shear_stress;
}

void mesh::save_boundary(const std::string& basename, bool verbose) {
    if (m_boundary_to_cell.empty() || m_offset_to_boundary.empty() || m_normals.empty()) extract_boundary();

    std::string boundary_name = basename + "_boundary.nrrd";
    std::string offset_name = basename + "_offset.nrrd";
    std::string normal_name = basename + "_normals.nrrd";
    std::vector<long int> array(5*m_boundary_to_cell.size());
    for (size_t i=0; i<m_boundary_to_cell.size(); i++) {
        const std::pair<face, index_type>& pfi = m_boundary_to_cell[i];
        const face& f = pfi.first;
        const index_type& n = pfi.second;
        array[5*i  ] = f[0];
        array[5*i+1] = f[1];
        array[5*i+2] = f[2];
        array[5*i+3] = f[3];
        array[5*i+4] = n;
    }
    size_t dims[2] = { 5, m_boundary_to_cell.size() };
    spurt::nrrd_utils::writeNrrd(&array[0], boundary_name, nrrdTypeLLong, 2, dims);
    if (verbose) {
        std::cout << "boundary to cell information written in "
        << boundary_name << "\n";
    }

    array.resize(m_offset_to_boundary.size()*5);
    for (size_t i=0; i<m_offset_to_boundary.size(); i++) {
        const std::pair<face, index_type>& pfi = m_offset_to_boundary[i];
        const face& f = pfi.first;
        const index_type& n = pfi.second;
        array[5*i  ] = f[0];
        array[5*i+1] = f[1];
        array[5*i+2] = f[2];
        array[5*i+3] = f[3];
        array[5*i+4] = n;
    }
    dims[1] = m_offset_to_boundary.size();
    spurt::nrrd_utils::writeNrrd(&array[0], offset_name, nrrdTypeLLong, 2, dims);
    if (verbose) {
        std::cout << "offset to boundary information written in "
        << offset_name << '\n';
    }

    dims[0] = 3;
    dims[1] = m_normals.size();
    spurt::nrrd_utils::writeNrrd(&m_normals[0], normal_name, nrrdTypeFloat, 2, dims);
    if (verbose) {
        std::cout << "boundary normals written in " << normal_name << '\n';
    }
}

void mesh::load_boundary(const std::string& basename, bool verbose) {
    std::string boundary_name = basename + "_boundary.nrrd";
    std::string offset_name = basename + "_offset.nrrd";
    std::string normal_name = basename + "_normals.nrrd";

    Nrrd* nin = spurt::nrrd_utils::readNrrd(boundary_name);
    m_boundary_to_cell.resize(nin->axis[1].size);
    long int* ids = (long int*)nin->data;
    for (size_t i=0; i<m_boundary_to_cell.size(); ++i) {
        face& f = m_boundary_to_cell[i].first;
        f[0] = ids[5*i  ];
        f[1] = ids[5*i+1];
        f[2] = ids[5*i+2];
        f[3] = ids[5*i+3];
        m_boundary_to_cell[i].second = ids[5*i+4];
    }
    if (verbose) {
        std::cout << "boundary to cell imported from "
        << boundary_name << "\n";
    }

    nin = spurt::nrrd_utils::readNrrd(offset_name);
    m_offset_to_boundary.resize(nin->axis[1].size);
    ids = (long int*)nin->data;
    for (size_t i=0; i<m_offset_to_boundary.size(); ++i) {
        face& f = m_offset_to_boundary[i].first;
        f[0] = ids[5*i  ];
        f[1] = ids[5*i+1];
        f[2] = ids[5*i+2];
        f[3] = ids[5*i+3];
        m_offset_to_boundary[i].second = ids[5*i+4];
    }
    if (verbose) {
        std::cout << "offset to boundary imported from "
        << offset_name << "\n";
    }

    nin = spurt::nrrd_utils::readNrrd(normal_name);
    m_normals.resize(nin->axis[1].size);
    float* nrmls = (float*)nin->data;
    for (size_t i=0; i<m_normals.size(); i++) {
        m_normals[i][0] = nrmls[3*i  ];
        m_normals[i][1] = nrmls[3*i+1];
        m_normals[i][2] = nrmls[3*i+2];
    }
    if (verbose) {
        std::cout << "boundary normals imported from "
        << normal_name << '\n';
    }
}

} // namespace spurt

#endif // __XAVIER_GEOMETRY_BOUNDARY_EXTRACTION_HPP__
