#ifndef __COHEN_SUTHERLAND_HPP__
#define __COHEN_SUTHERLAND_HPP__

#include <utility>
#include <math/types.hpp>
#include <math/bounding_box.hpp>

namespace cohen_sutherland {

typedef int     OutCode;

const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000

// Compute the bit code for a point (x, y) using the clip rectangle
// bounded diagonally by (xmin, ymin), and (xmax, ymax)
OutCode ComputeOutCode(double x, double y, const spurt::bbox2& bounds)
{
    OutCode code;
    
    const double& xmin = bounds.min()[0];
    const double& ymin = bounds.min()[1];
    const double& xmax = bounds.max()[0];
    const double& ymax = bounds.max()[1];
    
    code = INSIDE;          // initialised as being inside of clip window
    
    if (x < xmin) {         // to the left of clip window
        code |= LEFT;
    } else if (x > xmax) {  // to the right of clip window
        code |= RIGHT;
    }
    if (y < ymin) {         // below the clip window
        code |= BOTTOM;
    } else if (y > ymax) {  // above the clip window
        code |= TOP;
    }
    
    return code;
}

// Cohen–Sutherland clipping algorithm clips a line from
// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
// diagonal from (xmin, ymin) to (xmax, ymax).
bool clip(std::pair<spurt::vec2, spurt::vec2>& clipped,
          const spurt::vec2& p0, const spurt::vec2& p1,
          const spurt::bbox2& bounds)
{
    bool accept = false;
    
    double x0 = p0[0];
    double y0 = p0[1];
    double x1 = p1[0];
    double y1 = p1[1];
    const double& xmin = bounds.min()[0];
    const double& ymin = bounds.min()[1];
    const double& xmax = bounds.max()[0];
    const double& ymax = bounds.max()[1];
    
    // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
    OutCode outcode0 = ComputeOutCode(x0, y0, bounds);
    OutCode outcode1 = ComputeOutCode(x1, y1, bounds);
    
    while (true) {
        if (!(outcode0 | outcode1)) {     //logical or is 0. Trivially accept and get out of loop
            accept = true;
            break;
        } else if (outcode0 & outcode1) { //logical and is not 0. Trivially reject and get out of loop
            break;
        } else {
            // failed both tests, so calculate the line segment to clip
            // from an outside point to an intersection with clip edge
            double x, y;
            
            // At least one endpoint is outside the clip rectangle; pick it.
            OutCode outcodeOut = outcode0? outcode0 : outcode1;
            
            // Now find the intersection point;
            // use formulas y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
            if (outcodeOut & TOP) {           // point is above the clip rectangle
                x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
                y = ymax;
            } else if (outcodeOut & BOTTOM) { // point is below the clip rectangle
                x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
                y = ymin;
            } else if (outcodeOut & RIGHT) {  // point is to the right of clip rectangle
                y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
                x = xmax;
            } else if (outcodeOut & LEFT) {   // point is to the left of clip rectangle
                y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
                x = xmin;
            }
            // Now we move outside point to intersection point to clip
            // and get ready for next pass.
            if (outcodeOut == outcode0) {
                x0 = x;
                y0 = y;
                outcode0 = ComputeOutCode(x0, y0, bounds);
            } else {
                x1 = x;
                y1 = y;
                outcode1 = ComputeOutCode(x1, y1, bounds);
            }
        }
    }
    if (accept) {
        clipped.first = spurt::vec2(x0, y0);
        clipped.second = spurt::vec2(x1, y1);
    }
    
    return accept;
}

}


#endif