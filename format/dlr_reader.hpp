#ifndef __SPURT_DLR_READER_HPP_
#define __SPURT_DLR_READER_HPP_

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

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

class dlr_reader {

bool verbose;

inline void check_nc( int status ) {
    if( status != NC_NOERR )
        throw std::runtime_error(nc_strerror(status));
}

void load_vertices(std::vector<nvis::fvec3>& pos, int file) {
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

void load_cells(std::vector<long>& ids,
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

bool load_attributes(std::vector<double>& data, int file, const std::string& name) {
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
            //     if (strncmp("x_", __name, 2) == 0) {
            //         data.resize(nPoints*3);
            //         size_t start = 0, count = nPoints;
            //         ptrdiff_t stride = 1, imap = 3;
            //
            //         check_nc(nc_get_varm_double(file, i, &start, &count,
            //                                     &stride, &imap, &data[0]));
            //
            //         __name[0] = 'y';
            //         check_nc(nc_inq_varid(file, __name, &i));
            //
            //         check_nc(nc_get_varm_double(file, i, &start, &count,
            //                                     &stride, &imap, &data[1]));
            //
            //         __name[0] = 'z';
            //         check_nc(nc_inq_varid(file, __name, &i));
            //
            //         check_nc(nc_get_varm_double(file, i, &start, &count,
            //                                     &stride, &imap, &data[2]));
            //     }
            //     else {
                    data.resize(nPoints);
                    size_t start = 0, count = nPoints;
                    check_nc(nc_get_vara_double(file, i, &start, &count, &data[0]));
                // }
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

public:
    dlr_reader(const std::string& grid_file_name, const std::string& data_file_name)
            : _gfn(grid_file_name), _dfn(data_file_name), verbose(false) {}

    void read_mesh(bool boundary,
                   std::vector<nvis::fvec3>& vertices,
                   std::vector<long>& cell_indices,
                   std::vector<std::pair<cell_type, long> >& cell_types,
                   bool verbose=false) {
        this->verbose = verbose;
        int file, nPointsDimId;
        size_t nPoints;

        check_nc( nc_open( _gfn.c_str(), NC_NOWRITE, &file ) );
        check_nc( nc_inq_dimid( file, "no_of_points", &nPointsDimId ) );
        check_nc( nc_inq_dimlen( file, nPointsDimId, &nPoints ) );
        load_vertices(vertices, file);
        load_cells(cell_indices, cell_types, file, boundary);
        nc_close(file);
    }

    void read_data(const std::string& data_name, std::vector<double>& data,
                   bool verbose=false) {
        read_data_from_file(_dfn, data_name, data, verbose);
    }

    void read_data_from_file(const std::string& filename,
                             const std::string& data_name,
                             std::vector<double>& data, bool verbose=false) {
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
    void read_vector_data(const std::string& data_name, std::vector<Vec3>& data,
                          bool verbose=false) {
        read_vector_data_from_file<Vec3>(_dfn, data_name, data, verbose);
    }

    template<typename Vec3>
    void read_vector_data_from_file(const std::string& filename,
                                    const std::string& data_name,
                                    std::vector<Vec3>& data, bool verbose=false) {
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

private:
    std::string _gfn, _dfn;
};

}

#endif
