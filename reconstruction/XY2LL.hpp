#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <math/fixed_vector.hpp>
#include <util/timer.hpp>
#include <map>
#include <set>

namespace spurt {

nvis::vec2 simple_linear_regression(const std::vector<double>& x, const std::vector<double>& y)
{
    double xy_ = 0;
    double x_ = 0;
    double y_ = 0;
    double xsq_ = 0;
    for (int i=0 ; i<x.size() ; i++) {
        xy_ += x[i]*y[i];
        x_ += x[i];
        y_ += y[i];
        xsq_ += x[i]*x[i];
    }
    xy_ /= (double)x.size();
    x_ /= (double)x.size();
    y_ /= (double)x.size();
    xsq_ /= (double)x.size();
    double beta = (xy_ - x_*y_)/(xsq_ - x_*x_);
    double alpha = y_ - beta*x_;
    return nvis::vec2(alpha, beta);
}

class XY2LL {

    std::string read_string(std::fstream& in) {
        std::string str;
        char c;
        while (true) {
            in >> c;
            if (c == ',') {
                break;
            }
            str.append(1, c);
        }
        return str;
    }
    
    double phi(double r) const {
        return r*r*r;
    }
    
    
public:
    XY2LL(const std::string& filename, const nvis::ivec2& shift = nvis::ivec2(0,0) ) {
        std::fstream in(filename.c_str());
        
        std::string name, tmp;
        int x, y;
        double la, lo;
        char c;
        while (!in.eof()) {
            in >> c;
            if (c == '#') {
                std::getline(in, name);
                continue;
            }
            in.putback(c);
            name = read_string(in);
            names.push_back(name);
            in >> x >> c >> y >> c >> la >> c >> lo;
            coordXY.push_back(nvis::vec2(x+shift[0],y+shift[1]));
            coordLL.push_back(nvis::vec2(lo,la));
        }
        in.close();
        
        // Assume linear transformation between longitude and x (Mercator projection)
        std::vector<double> longs(coordLL.size()), xs(coordLL.size());
        for (int i=0 ; i<coordLL.size() ; ++i) {
            longs[i] = coordLL[i][0];
            xs[i] = coordXY[i][0];
        }
        nvis::vec2 lr = simple_linear_regression(longs, xs);
        m_long = lr[1];
        p_long = lr[0];
        
        // Assume following formula for Phi to Y: ln((1+sin(phi)/(1-sin(phi)))
        std::vector<double> merc(coordLL.size()), ys(coordLL.size());
        for (int i=0 ; i<coordLL.size() ; ++i) {
            double phi = coordLL[i][1];
            phi *= M_PI/180.;
            merc[i] = 0.5*log((1+sin(phi))/(1-sin(phi)));
            ys[i] = coordXY[i][1];
        }
        lr = simple_linear_regression(merc, ys);
        m_lat = lr[1];
        p_lat = lr[0];
    }
    
    nvis::vec2 toLoLa(const nvis::vec2& xy) const {
        double longitude = (xy[0]-p_long)/m_long;
        double tmp = (xy[1]-p_lat)/m_lat;
        double latitude = atan(sinh(tmp))*180./M_PI;
        return nvis::vec2(longitude, latitude);
    }
    
    nvis::vec2 toXY(const nvis::vec2& ll) const {
        double x = m_long*ll[0] + p_long;
        double phi = ll[1]*M_PI/180.;
        double y = m_lat*(0.5*log((1+sin(phi))/(1-sin(phi)))) + p_lat;
        return nvis::vec2(x, y);
    }
    
    const std::vector<std::string>& cities() const {
        return names;
    }
    
    double perimeter() const {
        return m_long*360.;
    }
    
    const std::vector<nvis::vec2>& xys() const {
        return coordXY;
    }
    
    const std::vector<nvis::vec2>& lls() const {
        return coordLL;
    }
    
private:
    std::vector<nvis::vec2> coordXY, coordLL;
    std::vector<std::string> names;
    double m_long, m_lat, p_long, p_lat;
};

};


