/*
 * Copyright 2002-2010 Guillaume Cottenceau.
 *
 * This software may be freely redistributed under the terms
 * of the X11 license.
 *
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <stdarg.h>
#include <data/raster.hpp>

#define PNG_DEBUG 3
#include <png.h>

namespace xavier {
    
class PNGWriter {
public:
    typedef nvis::fixed_vector<unsigned short, 3> color_t;
    typedef std::vector<color_t> row_t;
    
private:
    bool error_msg(const std::string& str) {
        std::cerr << "PNGWriter: ERROR: " << str << '\n';
        return false;
    }
    
    template<typename Raster_, size_t N>
    struct make_color {
        typedef typename Raster_::value_type     value_t;
        typedef typename Raster_::scalar_type    scalar_t;
        typedef typename Raster_::const_iterator const_iter_t;
        
        make_color(const Raster_ raster) {
            scalar_t minv=std::numeric_limits<scalar_t>::max();
            scalar_t maxv=std::numeric_limits<scalar_t>::min();
            for (const_iter_t it=raster.begin(); it!=raster.end(); ++it) {
                scalar_t _min = *std::min_element(it->begin(), it->end());
                scalar_t _max = *std::max_element(it->begin(), it->end());
                if (_min<minv) minv=_min;
                if (_max>maxv) maxv=_max;
            }
        }
        color_t operator()(const value_t& v) {
            unsigned short r=
                static_cast<unsigned short>((v[0]-minv)*65536/(maxv-minv));
            unsigned short g=
                static_cast<unsigned short>((v[1]-minv)*65536/(maxv-minv));
            if (N==3) {
                unsigned short b=
                    static_cast<unsigned short>(v[2]-minv*65536/(maxv-minv));
                return color_t(r, g, b);
            }
            else return color_t(r, 0, g);
        }
        scalar_t minv, maxv;
    };
    
    template<typename Raster_>
    struct make_color<Raster_, 1> {
        typedef typename Raster_::value_type  value_t;
        typedef typename Raster_::scalar_type scalar_t;
        
        make_color(const Raster_& raster) {
            minv=*std::min_element(raster.begin(), raster.end());
            maxv=*std::max_element(raster.begin(), raster.end());
        }
        color_t operator()(const value_t& v) {
            unsigned short u=
                static_cast<unsigned short>((v-minv)*65536/(maxv-minv));
            return color_t(u,u,u);
        }
        scalar_t minv, maxv;
    };
    
public:
    template<typename Raster_, size_t N=Raster_::dim>
    bool write(const Raster_& raster, const std::string& filename) {
        
        typedef typename Raster_::value_type     value_t;
        typedef typename Raster_::scalar_type    scalar_t;
        typedef typename Raster_::coord_type     coord_t;
        typedef typename Raster_::grid_type      grid_t;
        typedef typename Raster_::const_iterator const_iter_t;
        
        coord_t res=raster.grid().resolution();
        int width=res[0];
        int height=res[1];
        
        make_color<Raster_, N> colorize(raster);
        
        std::vector<row_t> image(height);
        for (int i=0; i<height; ++i) {
            image[i].resize(width);
            for (int j=0; j<width; ++j) {
                const value_t& v=raster(j,i);
                image[i][j]=colorize(v);
            }
        }
        
        FILE *fp = fopen(filename.c_str(), "wb");
        if (!fp) {
            return error_msg("WARNING: Unable to open " + filename + 
                             " for write operation");
        }
        
        // Initialization
        png_structp png_ptr = 
            png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!png_ptr) return error_msg("png_create_write_struct failed");
        
        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) return error_msg("png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr))) 
            return error_msg("Error during I/O initialization");

        png_init_io(png_ptr, fp);

        // Header
        if (setjmp(png_jmpbuf(png_ptr))) 
            return error_msg("Error during header export");

        png_set_IHDR(png_ptr, info_ptr, width, height,
                     16, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        png_write_info(png_ptr, info_ptr);

        // Write buffer to file
        if (setjmp(png_jmpbuf(png_ptr)))
            error_msg("Error during buffer export");
        
        png_set_swap(png_ptr);
        for (int i=0; i<height; ++i) {
            png_write_row(png_ptr, (png_bytep)(&image[i][0]));
        }
        png_write_end(png_ptr, NULL);
        fclose(fp);
        
        return true;
    }
};

} // namespace xavier
