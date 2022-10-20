#ifndef __FORMAT_DLRREADER_HPP__
#define __FORMAT_DLRREADER_HPP__

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include <math/fixed_vector.hpp>
#include <stdexcept>
#include <sstream>

#include <cassert>
#include <list>

#include <stdlib.h>
#include <string>

#include <fstream>
#include <unistd.h>
#include <netcdf.h>

namespace spurt {

class DLRreader {
    bool verbose;
    inline void check_nc( int status );
    void load_vertices(std::vector<nvis::fvec3>& pos, int file);
    struct dimLen {
        size_t nCells, nPointsPerCell;
        int varId;
    };

public:
    enum cell_type {
        NONE = 0,
        TRIANGLE,
        QUADRILATERAL,
        TETRAHEDRON,
        HEXAHEDRON,
        PRISM,
        PYRAMID,
    };

private:

    struct TCellConv {
        TCellConv(const char* _name, cell_type _type) : name(_name), type(_type) {}
        const char *const name;
        cell_type type;
    };

    void extract_boundary(std::vector<long>& newids,
                          std::vector<std::pair<cell_type, long> >& newoffsets,
                          const std::vector<long>& ids,
                          const std::vector<std::pair<cell_type, long> >& offsets);

    void load_cells(std::vector<long>& ids,
                    std::vector<std::pair<cell_type, long> >& offsets,
                    int file, bool boundary);

    bool load_attributes(std::vector<double>& data, int file, const std::string& name);

public:
    DLRreader(const std::string& grid_file_name, const std::string& data_file_name);

    void read_mesh(bool boundary,
                   std::vector<nvis::fvec3>& vertices,
                   std::vector<long>& cell_indices,
                   std::vector<std::pair<cell_type, long> >& cell_types,
                   bool verbose=false);

    void read_data(const std::string& data_name, std::vector<double>& data,
                   bool verbose=false);

    void read_data_from_file(const std::string& filename, const std::string& data_name,
                             std::vector<double>& data, bool verbose=false);

    template<typename Vec3>
    void read_vector_data(const std::string& data_name, std::vector<Vec3>& data,
                          bool verbose=false);

    template<typename Vec3>
    void read_vector_data_from_file(const std::string& filename, const std::string& data_name,
                                    std::vector<Vec3>& data, bool verbose=false);

private:
    std::string _gfn, _dfn;
};


inline void DLRreader::check_nc( int status ) {
    if( status != NC_NOERR )
        throw std::runtime_error(nc_strerror(status));
}

void DLRreader::load_vertices(std::vector<nvis::fvec3>& pos, int file) {
    int id, nPointsDimId;
    size_t nPoints;

    check_nc(nc_inq_dimid(file, "no_of_points", &nPointsDimId));
    check_nc(nc_inq_dimlen(file, nPointsDimId, &nPoints));

    std::vector<double> points(nPoints*3);

    size_t count = nPoints, start = 0;
    ptrdiff_t stride = 1, imap = 3;

    check_nc(nc_inq_varid(file, "points_xc", &id));
    check_nc(nc_get_varm_double(file, id, &start, &count,
                                &stride, &imap, &points[0]));

    check_nc(nc_inq_varid(file, "points_yc", &id));
    check_nc(nc_get_varm_double(file, id, &start, &count,
                                &stride, &imap, &points[1]));

    check_nc(nc_inq_varid(file, "points_zc", &id));
    check_nc(nc_get_varm_double(file, id, &start, &count,
                                &stride, &imap, &points[2]));

    if (verbose) std::cerr << "vertices successfully imported\n";

    pos.resize(nPoints);
    for (int i=0 ; i<nPoints ; ++i) {
        pos[i] = nvis::fvec3(points[3*i], points[3*i+1], points[3*i+2]);
    }
}

void DLRreader::extract_boundary(std::vector<long>& newids,
                                 std::vector<std::pair<cell_type, long> >& newoffsets,
                                 const std::vector<long>& ids,
                                 const std::vector<std::pair<cell_type, long> >& offsets) {

    struct face_index : public std::array<long, 4> {
        typedef std::array<long, 4> base_type;
        void sort(bool do_sort=true) {
            std::copy(this->begin(), this->end(), sorted.begin());
            if (do_sort)
                std::sort(sorted.begin(), sorted.end());
        }
        face_index() : base_type({-1, -1, -1, -1}) { sort(false); }
        face_index(long i, long j, long k, long l=-1)
            : base_type({-1, i, j, k}) { sort(); }
        face_index(const face_index& other) : base_type(other) { sort(); }
        std::array<long, 4> sorted;
    };

    struct Lt {
    bool operator()(const face_index& _this, const face_index& _that) const {
        for (int i=0; i<4; ++i) {
            if (_that.sorted[i] > _this.sorted[i]) return true;
            else if (_that.sorted[i] < _this.sorted[i]) return false;
        }
        return false;
    }};

    std::map<face_index, int, Lt> face_to_cell;
    std::vector<face_index> f(6);
    for (size_t i=0; i<offsets.size()-1; ++i)
    {
        cell_type cell = offsets[i].first;
        long start = offsets[i].second;
        long end = offsets[i+1].second;
        std::vector<long> v(ids.begin()+start, ids.begin()+end);
        int nfaces=0;
        if (cell == TETRAHEDRON) {
            f[0] = face_index(v[0], v[1], v[2]);
            f[1] = face_index(v[0], v[1], v[3]);
            f[2] = face_index(v[0], v[3], v[2]);
            f[3] = face_index(v[1], v[2], v[3]);
            nfaces = 4;
        }
        else if (cell == HEXAHEDRON) {
            f[0] = face_index(v[0], v[1], v[3], v[2]);
            f[1] = face_index(v[4], v[5], v[7], v[6]);
            f[2] = face_index(v[0], v[1], v[5], v[4]);
            f[3] = face_index(v[1], v[3], v[7], v[5]);
            f[4] = face_index(v[3], v[2], v[6], v[7]);
            f[5] = face_index(v[2], v[0], v[4], v[6]);
            nfaces = 6;
        }
        else if (cell == PRISM) {
            f[0] = face_index(v[0], v[1], v[2]);
            f[1] = face_index(v[4], v[3], v[5]);
            f[2] = face_index(v[0], v[2], v[5], v[3]);
            f[3] = face_index(v[1], v[4], v[5], v[2]);
            f[4] = face_index(v[0], v[1], v[4], v[3]);
            nfaces = 5;
        }
        else if (cell==PYRAMID) {
            f[0] = face_index(v[0], v[1], v[4]);
            f[1] = face_index(v[1], v[2], v[4]);
            f[2] = face_index(v[3], v[0], v[4]);
            f[3] = face_index(v[0], v[3], v[2], v[1]);
            nfaces = 4;
        }

        for (int j=0; j<nfaces; ++j) {
            auto iter = face_to_cell.find(f[j]);
            if (iter != face_to_cell.end()) {
                iter->second = iter->second + 1;
            }
            else face_to_cell[f[j]] = 1;
        }
    }

    if (verbose) {
        std::cout << "there are initially " << face_to_cell.size() << " unique faces in the entire mesh\n";
    }

    newids.clear();
    newoffsets.clear();
    size_t nskipped = 0;
    for (auto iter=face_to_cell.begin(); iter!=face_to_cell.end(); ++iter) {
        if (iter->second == 1) { // only one incident cell
            const face_index& fid = iter->first;
            if (fid.sorted[0] == -1) {
                if (fid.sorted[1] == -1) throw std::runtime_error("2 invalid indices in face");
                // TRIANGLE
                newoffsets.push_back(std::make_pair(TRIANGLE, newids.size())); // current offset
                newids.push_back(fid[1]);
                newids.push_back(fid[2]);
                newids.push_back(fid[3]);
            }
            else {
                // QUADRILATERAL
                newoffsets.push_back(std::make_pair(QUADRILATERAL, newids.size())); // current offset
                newids.push_back(fid[0]);
                newids.push_back(fid[1]);
                newids.push_back(fid[2]);
                newids.push_back(fid[3]);
            }
        }
        else ++nskipped;
    }
    newoffsets.push_back(std::make_pair(NONE, newids.size())); // add an offset past the last index

    if (verbose) {
        std::cout << newoffsets.size()-1 << " boundary cells found\n";
        std::cout << nskipped << " internal faces\n";
    }
}

void DLRreader::load_cells(std::vector<long>& ids,
                std::vector<std::pair<cell_type, long> >& offsets,
                int file, bool boundary) {
    int allCells = 0, allCellsSize = 0, nPointsDimId;
    size_t nPoints;

    static const TCellConv cellCore[5] = {
        TCellConv("points_of_hexaeders", HEXAHEDRON),
        TCellConv("points_of_tetraeders", TETRAHEDRON),
        TCellConv("points_of_prisms", PRISM),
        TCellConv("points_of_pyramids", PYRAMID),
        TCellConv(0, NONE) //delimiter
    };

    static const TCellConv cellBdry[3] = {
        TCellConv("points_of_surfacetriangles", TRIANGLE),
        TCellConv("points_of_surfacequadrilaterals", QUADRILATERAL),
        TCellConv(0, NONE) //delimiter
    };

    std::vector<dimLen> dimLens;

    const TCellConv *cellConv = boundary ? cellBdry : cellCore;

    check_nc(nc_inq_dimid(file, "no_of_points", &nPointsDimId));
    check_nc(nc_inq_dimlen(file, nPointsDimId, &nPoints));

    for (int i = 0; cellConv[i].name != 0; i++) {
        dimLen dl;
        dl.nCells = 0;

        if (nc_inq_varid(file, cellConv[i].name,
                         &dl.varId) == NC_NOERR) {
            int dimIds[2];

            check_nc(nc_inq_vardimid(file, dl.varId, dimIds));
            check_nc(nc_inq_dimlen(file, dimIds[0], &dl.nCells));
            check_nc(nc_inq_dimlen(file, dimIds[1], &dl.nPointsPerCell));

            allCells += dl.nCells;
            allCellsSize += dl.nPointsPerCell * dl.nCells;
        }
        dimLens.push_back(dl);
    }

    ids.resize(allCellsSize);
    offsets.resize(allCells+1);

    typedef std::vector<std::pair<cell_type, long> > typeVector;

    typeVector::iterator k = offsets.begin(), kend;

    unsigned int actInd = 0;

    for (int i = 0; cellConv[i].name != 0; i++) {
        if (dimLens[i].nCells) {
            size_t start[] = { 0, 0 };
            size_t count[] = { dimLens[i].nCells, dimLens[i].nPointsPerCell };

            char name[256];
            nc_inq_varname(file, dimLens[i].varId, name);

            check_nc(nc_get_vara_long(file, dimLens[i].varId,
                                      start, count,
                                      (long*)&ids[actInd]));

            cell_type t = cellConv[i].type;
            kend = k + dimLens[i].nCells;

            for (; k != kend; k++, actInd += dimLens[i].nPointsPerCell) {
                k->first = t;
                k->second = actInd;
            }
        }
    }

    k->first = NONE;
    k->second = actInd;

    if (verbose) std::cerr << "cell definitions successfully imported\n";
}

bool DLRreader::load_attributes(std::vector<double>& data, int file, const std::string& name) {
    int nvars, nPointsDimId;
    size_t nPoints;

    check_nc(nc_inq_nvars(file, &nvars));
    check_nc(nc_inq_dimid(file, "no_of_points", &nPointsDimId));
    check_nc(nc_inq_dimlen(file, nPointsDimId, &nPoints));

    if (verbose) std::cerr << "there are " << nvars << " fields in input.\n";
    if (verbose) std::cout << "requested field: " << name << '\n';

    for (int i = 0; i < nvars; i++) {
        char __name[256];

        check_nc(nc_inq_varname(file, i, __name));

        if (strcmp(__name, name.c_str())) {
            if (verbose) std::cout << "skipping attribute " << __name << '\n';
            continue;
        }
        int ndims;

        check_nc(nc_inq_varndims(file, i, &ndims));

        assert(ndims < 1000);
        assert(ndims > 0);

        int *dimid = new int[ndims];

        check_nc(nc_inq_vardimid(file, i, &(dimid[0])));

        if (dimid[0] == nPointsDimId) {
            if (ndims == 1) {
                data.resize(nPoints);
                size_t start = 0, count = nPoints;
                check_nc(nc_get_vara_double(file, i, &start, &count, &data[0]));
            }
            else {
                size_t tsize = 1, tdim;
                size_t *start = new size_t[ndims];
                size_t *count = new size_t[ndims];

                check_nc(nc_inq_dimlen(file, dimid[1], &tdim));

                for (int j = 1; j < ndims; ++j) {
                    size_t tmp;
                    check_nc(nc_inq_dimlen(file, dimid[j], &tmp));
                    assert(tmp == tdim);
                    tsize *= tdim;

                    start[j] = 0;
                    count[j] = tdim;
                }

                start[0] = 0;
                count[0] = nPoints;

                data.resize(nPoints*tsize);

                check_nc(nc_get_vara_double(file, i, &start[0], &count[0] , &data[0]));
                delete [] start;
                delete [] count;
            }
        }
        delete [] dimid;
        return true;
    }
    return false;
}

DLRreader::DLRreader(const std::string& grid_file_name, const std::string& data_file_name)
            : _gfn(grid_file_name), _dfn(data_file_name), verbose(false) {}

void DLRreader::read_mesh(bool boundary,
               std::vector<nvis::fvec3>& vertices,
               std::vector<long>& cell_indices,
               std::vector<std::pair<cell_type, long> >& cell_types,
               bool verbose) {
    this->verbose = verbose;
    int file, nPointsDimId;
    size_t nPoints;

    check_nc( nc_open( _gfn.c_str(), NC_NOWRITE, &file ) );
    check_nc( nc_inq_dimid( file, "no_of_points", &nPointsDimId ) );
    check_nc( nc_inq_dimlen( file, nPointsDimId, &nPoints ) );
    load_vertices(vertices, file);
    std::cout << vertices.size() << " vertices imported\n";

    load_cells(cell_indices, cell_types, file, boundary);
    if (boundary && cell_indices.empty()) {
        if (verbose) {
            std::cout << "requested boundary was not found in file. Extracting...\n";
        }
        std::vector<long> tmp_ids;
        std::vector<std::pair<cell_type, long> > tmp_types;
        // no boundary mesh in DataSet
        load_cells(tmp_ids, tmp_types, file, false);
        extract_boundary(cell_indices, cell_types, tmp_ids, tmp_types);
    }
    nc_close(file);
}

void DLRreader::read_data(const std::string& data_name, std::vector<double>& data,
                          bool verbose) {

    read_data_from_file(_dfn, data_name, data, verbose);
}

void DLRreader::read_data_from_file(const std::string& filename,
                        const std::string& data_name,
                         std::vector<double>& data, bool verbose) {
    this->verbose = verbose;
    int file, nPointsDimId;
    size_t nPoints;

    check_nc( nc_open(filename.c_str(), NC_NOWRITE, &file ) );
    if (verbose) std::cout << "file=" << file << std::endl;
    if (!load_attributes(data, file, data_name)) {
        std::cerr << "requested data is not available" << std::endl;
    }
    nc_close(file);
}

template<typename Vec3>
void DLRreader::read_vector_data(const std::string& data_name, std::vector<Vec3>& data,
                                 bool verbose) {
    read_vector_data_from_file<Vec3>(_dfn, data_name, data, verbose);
}

template<typename Vec3>
void DLRreader::read_vector_data_from_file(const std::string& filename, const std::string& data_name,
                                           std::vector<Vec3>& data, bool verbose) {
    this->verbose = verbose;
    int file, nPointsDimId;
    size_t nPoints;

    check_nc( nc_open(filename.c_str(), NC_NOWRITE, &file ) );
    check_nc( nc_inq_dimid( file, "no_of_points", &nPointsDimId ) );
    check_nc( nc_inq_dimlen( file, nPointsDimId, &nPoints ) );
    std::vector<double> x, y, z;
    if (!load_attributes(x, file, "x_" + data_name)) {
        std::cerr << "requested data (x_" << data_name << ") is not available" << std::endl;
    }
    if (!load_attributes(y, file, "y_" + data_name)) {
        std::cerr << "requested data (y_" << data_name << ") is not available" << std::endl;
    }
    if (!load_attributes(z, file, "z_" + data_name)) {
        std::cerr << "requested data (z_" << data_name << ") is not available" << std::endl;
    }
    nc_close(file);
    data.resize(x.size());
    for (size_t i=0; i<x.size(); ++i) {
        data[i][0] = x[i];
        data[i][1] = y[i];
        data[i][2] = z[i];
    }
}

} // namespace spurt

#endif // __FORMAT_DLRREADER_HPP__
