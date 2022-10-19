#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <math.h>

// to do:  add slerp

namespace MathUtils
{

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

#ifndef M_2_PI
#define M_2_PI 6.28318530718
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679
#endif

#ifndef M_PI_4
#define M_PI_4 0.785398163397
#endif

#ifndef M_PI_8
#define M_PI_8 0.392699081699
#endif

#ifndef M_e
#define M_e 2.71828182846
#endif

#ifndef RADIANS_TO_DEGREES
#define RADIANS_TO_DEGREES 180.0/M_PI
#endif

#ifndef DEGREES_TO_RADIANS
#define DEGREES_TO_RADIANS M_PI/180.0
#endif

static inline double radiansToDegrees(const double radians)
{
	return radians*RADIANS_TO_DEGREES;
}
static inline double degreesToRadians(const double degrees)
{
	return degrees*DEGREES_TO_RADIANS;
}

template <class T>
static inline T sqr(const T &a)
{
	return a*a;
}

template <class T> static inline T lerp(const double t, const T x1, const T x2)
{
	return (1.0 - t)*x1 + t*x2;
}
template<> inline unsigned char lerp(const double t, const unsigned char x1, const unsigned char x2)
{
	return (unsigned char)((1.0 - t)*(float)x1 + t*(float)x2);
}

template <class T> static inline T bilinearLerp(const T val00, const T val01, const T val10, const T val11, const double tx, const double ty)
{
	typename T::value_type u(tx), v(ty);

	T i1 = v * val00 + (1 - v) * val01;
	T i2 = v * val10 + (1 - v) * val11;
	return u*i1 + (1 - u)*i2;
}
template <> inline unsigned char bilinearLerp(const unsigned char val00, const unsigned char val01, const unsigned char val10,
        const unsigned char val11, const double tx, const double ty)
{
	double ty2 = 1.0 - ty;
	double i1 = ty * (double)val00 + ty2 * (double)val01;
	double i2 = ty * (double)val10 + ty2 * (double)val11;
	return (unsigned char)(tx*i1 + (1.0 - tx)*i2);
}

template <class T> static inline T bilinearLerp(const T val00, const T val01, const T val10, const T val11,
        const double xUp,  const double yUp, const double xMid, const double yMid)
{
	return bilinearLerp(val00, val01, val10, val11, xUp - xMid, yUp - yMid);
}

template <class T> static inline T trilinearLerp(const T val000, const T val001, const T val010, const T val011,
        const T val100, const T val101, const T val110, const T val111,
        const double tx, const double ty, const double tz)
{
	double tz2 = 1.0 - tz;
	double ty2 = 1.0 - ty;

	T i1 = tz * val000 + tz2 * val001;
	T i2 = tz * val100 + tz2 * val101;
	T j1 = tz * val010 + tz2 * val011;
	T j2 = tz * val110 + tz2 * val111;

	T w1 = ty * i1 + ty2 * j1;
	T w2 = ty * i2 + ty2 * j2;

	return tx*w1 + (1.0 - tx)*w2;
}
template <> inline unsigned char trilinearLerp(const unsigned char val000, const unsigned char val001, const unsigned char val010, const unsigned char val011,
        const unsigned char val100, const unsigned char val101, const unsigned char val110, const unsigned char val111,
        const double tx, const double ty, const double tz)
{
	double tz2 = 1.0 - tz;
	double ty2 = 1.0 - ty;

	double i1 = tz * (double)val000 + tz2 * (double)val001;
	double i2 = tz * (double)val100 + tz2 * (double)val101;
	double j1 = tz * (double)val010 + tz2 * (double)val011;
	double j2 = tz * (double)val110 + tz2 * (double)val111;

	double w1 = ty * i1 + ty2 * j1;
	double w2 = ty * i2 + ty2 * j2;

	return (unsigned char)(tx*w1 + (1.0 - tx)*w2);
}

template <class T> static inline T trilinearLerp(const T val000, const T val001, const T val010, const T val011,
        const T val100, const T val101, const T val110, const T val111,
        const double xUp,  const double yUp,  const double zUp,
        const double xMid, const double yMid, const double zMid)
{
	return trilinearLerp(val000, val001, val010, val011, val100, val101, val110, val111, xUp - xMid, yUp - yMid, zUp - zMid);
}

// 1D interpolation using Catmull-Rom spline
// t should be [0, 1] and describes a point between q2 and q3, where t = 0 -> q2  and t = 1 -> q3
template <class T> static inline T catmullRomInterpolate(const T q1, const T q2, const T q3, const T q4, const double t,
        const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	// q1 = q(i-1)
	// q2 = q(i)
	// q3 = q(i+1)
	// q4 = q(i+2)
	T d1 = (q3 - q1) * 0.5;
	T d2 = (q4 - q2) * 0.5;
	T dq = q3 - q2;

	// if q(i) < q(i+1) and d2 < 0, the max of cubic occurs somewhere in the middle and is greater than q(i+1)
	//  this can lead to creating a return value that is greater than any of the q's which can lead to
	//  expontential increases - by checking the signs we can prevent this
	if (preserveMonotonicity && !clampToLocalBounds) {
		// one fix is to set the slopes to zero if they have different signs than delta(q)
		//  but this smothing may be a bit excessive
		if (dq < 0) {
			if (d1 > 0) d1 = 0;
			if (d2 > 0) d2 = 0;
		}
		else {
			if (d1 < 0) d1 = 0;
			if (d2 < 0) d2 = 0;
		}
	}

	T result = q2 + t * (d1 + t * (3.0 * dq - 2.0 * d1 - d2 + t * (-2.0 * dq + d1 + d2)));

	// another fix is to clamp against q(i) and q(i+1), but can result in sharp discontinuities
	if (preserveMonotonicity && clampToLocalBounds) {
		if (q2 < q3) result = clamp(result, q2, q3);
		else result = clamp(result, q3, q2);
	}

	return result;
}

// 1D interpolation using Catmull-Rom spline - computes t from provided values (assumes computed t will be [0, 1])
template <class T> static inline T catmullRomInterpolate(const T q1, const T q2, const T q3, const T q4, double const q2Loc, const double interpLoc,
        const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	return catmullRomInterpolate(q1, q2, q3, q4, interpLoc - q2Loc, preserveMonotonicity, clampToLocalBounds);
}

// 1D interpolation using Catmull-Rom spline - computes t from provided values (normalizes t to [0, 1])
template <class T> static inline T catmullRomInterpolate(const T q1, const T q2, const T q3, const T q4,
        const double q2Loc, const double interpLoc, const double q3Loc,
        const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	return catmullRomInterpolate(q1, q2, q3, q4, (interpLoc - q2Loc) / (q3Loc - q2Loc), preserveMonotonicity, clampToLocalBounds);
}

// 2D interpolation using Catmull-Rom spline
// t should be [0, 1] and describes a point between q2 and q3, where t = 0 -> q2  and t = 1 -> q3
template <class T> static inline T bicubicCatmullRomInterpolate(const T q11, const T q12, const T q13, const T q14, const T q21, const T q22, const T q23, const T q24,
        const T q31, const T q32, const T q33, const T q34, const T q41, const T q42, const T q43, const T q44,
        const double tx, const double ty, const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	return catmullRomInterpolate(catmullRomInterpolate(q11, q12, q13, q14, ty, preserveMonotonicity, clampToLocalBounds),
	                             catmullRomInterpolate(q21, q22, q23, q24, ty, preserveMonotonicity, clampToLocalBounds),
	                             catmullRomInterpolate(q31, q32, q33, q34, ty, preserveMonotonicity, clampToLocalBounds),
	                             catmullRomInterpolate(q41, q42, q43, q44, ty, preserveMonotonicity, clampToLocalBounds), tx);
}

// 2D interpolation using Catmull-Rom spline, assumes q3Loc > interpLoc > q2Loc and q3Loc-q2Loc = 1
template <class T> static inline T bicubicCatmullRomInterpolate(const T q11, const T q12, const T q13, const T q14, const T q21, const T q22, const T q23, const T q24,
        const T q31, const T q32, const T q33, const T q34, const T q41, const T q42, const T q43, const T q44,
        const double q2LocX, const double q2LocY, const double interpLocX, const double interpLocY,
        const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	return tricubicCatmullRomInterpolate(q11, q12, q13, q14, q21, q22, q23, q24, q31, q32, q33, q34, q41, q42, q43, q44,
	                                     interpLocX - q2LocX, interpLocY - q2LocY, preserveMonotonicity, clampToLocalBounds);
}

// 3D interpolation using Catmull-Rom spline
// t should be [0, 1] and describes a point between q2 and q3, where t = 0 -> q2  and t = 1 -> q3
template <class T> static inline T tricubicCatmullRomInterpolate(const T q111, const T q112, const T q113, const T q114, const T q121, const T q122, const T q123, const T q124,
        const T q131, const T q132, const T q133, const T q134, const T q141, const T q142, const T q143, const T q144,
        const T q211, const T q212, const T q213, const T q214, const T q221, const T q222, const T q223, const T q224,
        const T q231, const T q232, const T q233, const T q234, const T q241, const T q242, const T q243, const T q244,
        const T q311, const T q312, const T q313, const T q314, const T q321, const T q322, const T q323, const T q324,
        const T q331, const T q332, const T q333, const T q334, const T q341, const T q342, const T q343, const T q344,
        const T q411, const T q412, const T q413, const T q414, const T q421, const T q422, const T q423, const T q424,
        const T q431, const T q432, const T q433, const T q434, const T q441, const T q442, const T q443, const T q444,
        const double tx, const double ty, const double tz,
        const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	return catmullRomInterpolate(
	           catmullRomInterpolate(catmullRomInterpolate(q111, q112, q113, q114, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q121, q122, q123, q124, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q131, q132, q133, q134, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q141, q142, q143, q144, tz, preserveMonotonicity, clampToLocalBounds), ty),
	           catmullRomInterpolate(catmullRomInterpolate(q211, q212, q213, q214, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q221, q222, q223, q224, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q231, q232, q233, q234, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q241, q242, q243, q244, tz, preserveMonotonicity, clampToLocalBounds), ty),
	           catmullRomInterpolate(catmullRomInterpolate(q311, q312, q313, q314, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q321, q322, q323, q324, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q331, q332, q333, q334, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q341, q342, q343, q344, tz, preserveMonotonicity, clampToLocalBounds), ty),
	           catmullRomInterpolate(catmullRomInterpolate(q411, q412, q413, q414, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q421, q422, q423, q424, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q431, q432, q433, q434, tz, preserveMonotonicity, clampToLocalBounds),
	                                 catmullRomInterpolate(q441, q442, q443, q444, tz, preserveMonotonicity, clampToLocalBounds), ty), tx);
}

// 3D interpolation using Catmull-Rom spline, assumes q3Loc > interpLoc > q2Loc and q3Loc-q2Loc = 1
template <class T> static inline T tricubicCatmullRomInterpolate(const T q111, const T q112, const T q113, const T q114, const T q121, const T q122, const T q123, const T q124,
        const T q131, const T q132, const T q133, const T q134, const T q141, const T q142, const T q143, const T q144,
        const T q211, const T q212, const T q213, const T q214, const T q221, const T q222, const T q223, const T q224,
        const T q231, const T q232, const T q233, const T q234, const T q241, const T q242, const T q243, const T q244,
        const T q311, const T q312, const T q313, const T q314, const T q321, const T q322, const T q323, const T q324,
        const T q331, const T q332, const T q333, const T q334, const T q341, const T q342, const T q343, const T q344,
        const T q411, const T q412, const T q413, const T q414, const T q421, const T q422, const T q423, const T q424,
        const T q431, const T q432, const T q433, const T q434, const T q441, const T q442, const T q443, const T q444,
        const double q2LocX, const double q2LocY, const double q2LocZ,
        const double interpLocX, const double interpLocY, const double interpLocZ,
        const bool preserveMonotonicity = false, const bool clampToLocalBounds = false)
{
	return tricubicCatmullRomInterpolate(q111, q112, q113, q114, q121, q122, q123, q124, q131, q132, q133, q134, q141, q142, q143, q144,
	                                     q211, q212, q213, q214, q221, q222, q223, q224, q231, q232, q233, q234, q241, q242, q243, q244,
	                                     q311, q312, q313, q314, q321, q322, q323, q324, q331, q332, q333, q334, q341, q342, q343, q344,
	                                     q411, q412, q413, q414, q421, q422, q423, q424, q431, q432, q433, q434, q441, q442, q443, q444,
	                                     interpLocX - q2LocX, interpLocY - q2LocY, interpLocZ - q2LocZ, preserveMonotonicity, clampToLocalBounds);
}

#ifndef round
#define round(a) ((int)(a+0.4999999999))
#endif

//#ifndef ceil
//	#define ceil(a) ((int)(a+0.9999999999))
//#endif

//#ifndef floor
//	#define floor(a) ((int)a)
//#endif

template <class T>
static inline T clamp(const T &val, const T &low, const T &high)
{
	if (val <= low) return low;
	else if (val >= high) return high;
	else return val;
}

template <class T>	// shortcuts clamping 2 vals with the relation val1 <= val2
static inline void clamp2(T &val1, T &val2, const T &low, const T &high)
{
	if (val2 <= low) val1 = val2 = low;
	else if (val1 >= high) val1 = val2 = high;
	else {
		if (val1 < low) val1 = low;
		if (val2 > high) val2 = high;
	}
}

static inline int mod(const int a, const int b)    // supports mod (%) for negative numbers
{
	int n = a / b;
	int c = a - n * b;
	if (c < 0) c += b;
	return c;
}

static inline void getBoundingIntegers(const double val, int &lower, int &upper, const double tolerance = 1e-8)
{
	lower = (int)val;
	if (val < 0) lower--;
	if (fabs(lower - val) < tolerance) upper = lower; // val = x.00001 or -x.99999
	else {
		upper = lower + 1;
		if (fabs(upper - val) < tolerance) lower = upper; // val = x.99999 or -x.00001
	}
}

// to do: make max/min static inline and in such a way that they don't cause compile errors with std::max
//          also so they don't have to be named max3, max4, etc.

//#ifndef max
//	#define max(a,b)(((a)>(b)) ? (a) : (b))
//#endif

#ifndef max3
template <class T>
static inline T max3(const T &a, const T &b, const T &c)
{
	return max(a, max(b, c));
}
#endif
template <class T>
static inline T max4(const T &a, const T &b, const T &c, const T &d)
{
	return max(max(a, b), max(c, d));
}

//#ifndef min
//	#define min(a,b)(((a)<(b)) ? (a) : (b))
//#endif

#ifndef min3
template <class T>
static inline T min3(const T &a, const T &b, const T &c)
{
	return min(a, min(b, c));
}
#endif
template <class T>
static inline T min4(const T &a, const T &b, const T &c, const T &d)
{
	return min(min(a, b), min(c, d));
}

template<class T>
inline void minMax(const T &a, T b, T &minVal, T &maxVal)
{
	if (a < b) {
		minVal = a;
		maxVal = b;
	}
	else {
		minVal = b;
		maxVal = a;
	}
}

template<class T>
inline void minMax(const T &a, const T &b, const T &c, T &minVal, T &maxVal)
{
	if (a < b) {
		if (a < c) {
			minVal = a;
			maxVal = max(b, c);
		}
		else {
			minVal = c;
			maxVal = max(a, b);
		}
	}
	else {
		if (b < c) {
			minVal = b;
			maxVal = max(a, c);
		}
		else {
			minVal = c;
			maxVal = a;
		}
	}
}

template<class T>
inline void minMax(const T &a, const T &b, const T &c, const T &d, T &minVal, T &maxVal)
{
	if (a < b) {
		if (c < d) {
			minVal = min(a, c);
			maxVal = max(b, d);
		}
		else {
			minVal = min(a, d);
			maxVal = max(b, c);
		}
	}
	else {
		if (c < d) {
			minVal = min(b, c);
			maxVal = max(a, d);
		}
		else {
			minVal = min(b, d);
			maxVal = max(a, c);
		}
	}
}

template<class T>
inline T smoothStep(const T &r)
{
	if (r < 0) return 0;
	else if (r > 1) return 1;
	return r*r*r*(10 + r*(-15 + r*6));
}

template<class T>
inline T smoothStep(const T &r, const T &rLower, const T &rUpper, const T &valLower, const T &valUpper)
{
	return valLower + smoothStep((r - rLower) / (rUpper - rLower)) * (valUpper - valLower);
}

template<class T>
inline T ramp(const T &r)
{
	return smoothStep((r + 1) / 2)*2 - 1;
}

static inline double getLog(double x, double logBase)
{
	return log(x) / log(logBase);
}
static inline double log2(const double x)
{
	return getLog(x, 2);
}
static inline int log2Int(const float x)
{
	return ((*(int*) &x) >> 23) - 127;
}

static inline bool isPowerOf2(const int x)
{
	return (x & (x - 1)) == 0;
}

static inline unsigned int roundUpPowerOf2(const unsigned int x)
{
	unsigned int y = x - 1;
	y |= y >> 1;
	y |= y >> 2;
	y |= y >> 4;
	y |= y >> 8;
	y |= y >> 16;
	return y + 1;
}

static inline unsigned int roundDownPowerOf2(const unsigned int x)
{
	unsigned int y = x;
	int exponent = 0;
	while (y > 1) {
		++exponent;
		y >>= 1;
	}
	return 1 << exponent;
}

// most significant bit = integer of log base 2 of x
static inline unsigned int getMostSignificantBit(const unsigned int x)   // O(lg n)
{
	unsigned int val = x;

	unsigned int result = (val > 0xFFFF) << 4;
	val >>= result;
	unsigned int shift  = (val > 0xFF) << 3;
	val >>= shift;
	result |= shift;
	shift  = (val > 0xF) << 2;
	val >>= shift;
	result |= shift;
	shift  = (val > 0x3) << 1;
	val >>= shift;
	result |= shift;
	return result | (val >> 1);
}

// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
static inline unsigned int getNumberOfTrailingZeros_Mod(const int x)
{
	static const int Mod37BitPosition[] = { // map a bit value mod 37 to its position	32, 0, 1, 26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11, 0, 13, 4, 7, 17,
		0, 25, 22, 31, 15, 29, 10, 12, 6, 0, 21, 14, 9, 5, 20, 8, 19, 18
	};
	return Mod37BitPosition[(-x & x) % 37];
}
static inline unsigned int getNumberOfTrailingZeros_Mult(const int x)
{
	static const int MultiplyDeBruijnBitPosition[32] = {	0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
	        31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
	                                                   };
	return MultiplyDeBruijnBitPosition[((x & -x) * 0x077CB531U) >> 27];
}

template<class T> static inline int getBit(const T x, const int bitIdx)
{
	if (x & 1 << bitIdx) return 1;
	else return 0;
}
template<class T> static inline void printBit(const T x, const int bitIdx)
{
	printf("%d", getBit(x, bitIdx));
}
template<class T> static inline void printBinary(const T x)
{
	//for (unsigned int mask = 0x80000000; mask != 0; mask >>= 1) printf("%c", (mask & x) ? '1' : '0');
	// to do: below is untested - used with the template, above works for int
	unsigned int numBits = sizeof(x) * 8;
	for (int i = numBits - 1; i >= 0; i--) printBit(x, i);
}
template<class T> static inline void zeroBit(T &x, const int bitIdx)
{
	x &= ~(1 << bitIdx);
}
template<class T> static inline void setBit(T &x, const int bitIdx)
{
	x |= 1 << bitIdx;
}

// to do: make more efficient, perhaps using a lookup table
// (Graphics Gems 1 has a chapter about this, "Bit Interleaving for Quad- or Octrees")
static inline unsigned int interleaveBits(const unsigned int x, const unsigned int y)
{
	// 16 shifts, 32 or's, 32 and's
	unsigned int val = y & 32768;
	val = (val << 1) | (x & 32768) | (y & 16384);
	val = (val << 1) | (x & 16384) | (y & 8192);
	val = (val << 1) | (x & 8192)  | (y & 4096);
	val = (val << 1) | (x & 4096)  | (y & 2048);
	val = (val << 1) | (x & 2048)  | (y & 1024);
	val = (val << 1) | (x & 1024)  | (y & 512);
	val = (val << 1) | (x & 512)   | (y & 256);
	val = (val << 1) | (x & 256)   | (y & 128);
	val = (val << 1) | (x & 128)   | (y & 64);
	val = (val << 1) | (x & 64)    | (y & 32);
	val = (val << 1) | (x & 32)    | (y & 16);
	val = (val << 1) | (x & 16)    | (y & 8);
	val = (val << 1) | (x & 8)     | (y & 4);
	val = (val << 1) | (x & 4)     | (y & 2);
	val = (val << 1) | (x & 2)     | (y & 1);
	return (val << 1) | (x&1);

	// 31 shifts, 32 or's, 32 and's
	/*unsigned int val = (x&1) | ((y&1) << 1);
	val |= ((x&2) << 1) | ((y&2) << 2);
	val |= ((x&4) << 2) | ((y&4) << 3);
	val |= ((x&8) << 3) | ((y&8) << 4);
	val |= ((x&16) << 4) | ((y&16) << 5);
	val |= ((x&32) << 5) | ((y&32) << 6);
	val |= ((x&64) << 6) | ((y&64) << 7);
	val |= ((x&128) << 7) | ((y&128) << 8);
	val |= ((x&256) << 8) | ((y&256) << 9);
	val |= ((x&512) << 9) | ((y&512) << 10);
	val |= ((x&1024) << 10) | ((y&1024) << 11);
	val |= ((x&2048) << 11) | ((y&2048) << 12);
	val |= ((x&4096) << 12) | ((y&4096) << 13);
	val |= ((x&8192) << 13) | ((y&8192) << 14);
	val |= ((x&16384) << 14) | ((y&16384) << 15);
	val |= ((x&32768) << 15) | ((y&32768) << 16);
	return val;*/
}
static inline void unInterleaveBits(const unsigned int val, unsigned int &x, unsigned int &y)
{
	x = (val & 0x1)              | ((val & 0x4)       >> 1)  | ((val & 0x10)       >> 2)  | ((val & 0x40)       >> 3)  |
	    ((val & 0x100)     >> 4)  | ((val & 0x400)     >> 5)  | ((val & 0x1000)     >> 6)  | ((val & 0x4000)     >> 7)  |
	    ((val & 0x10000)   >> 8)  | ((val & 0x40000)   >> 9)  | ((val & 0x100000)   >> 10) | ((val & 0x400000)   >> 11) |
	    ((val & 0x1000000) >> 12) | ((val & 0x4000000) >> 13) | ((val & 0x10000000) >> 14) | ((val & 0x40000000) >> 15);

	y = ((val & 0x2)       >> 1)  | ((val & 0x8)       >> 2)  | ((val & 0x20)       >> 3)  | ((val & 0x80)       >> 4)  |
	    ((val & 0x200)     >> 5)  | ((val & 0x800)     >> 6)  | ((val & 0x2000)     >> 7)  | ((val & 0x8000)     >> 8)  |
	    ((val & 0x20000)   >> 9)  | ((val & 0x80000)   >> 10) | ((val & 0x200000)   >> 11) | ((val & 0x800000)   >> 12) |
	    ((val & 0x2000000) >> 13) | ((val & 0x8000000) >> 14) | ((val & 0x20000000) >> 15) | ((val & 0x80000000) >> 16);
}

static inline unsigned int interleaveBits(const unsigned int x, const unsigned int y, const unsigned int z)
{
	// 20 shifts, 29 or's, 30 and's
	unsigned int val = z & 512;
	val = (val << 2) | ((z & 256) | (x & 512)) | ((y & 512) << 1);
	val = (val << 2) | ((z & 128) | (x & 256)) | ((y & 256) << 1);
	val = (val << 2) | ((z & 64)  | (x & 128)) | ((y & 128) << 1);
	val = (val << 2) | ((z & 32)  | (x & 64)) | ((y & 64)  << 1);
	val = (val << 2) | ((z & 16)  | (x & 32)) | ((y & 32)  << 1);
	val = (val << 2) | ((z & 8)   | (x & 16)) | ((y & 16)  << 1);
	val = (val << 2) | ((z & 4)   | (x & 8)) | ((y & 8)   << 1);
	val = (val << 2) | ((z & 2)   | (x & 4)) | ((y & 4)   << 1);
	val = (val << 2) | ((z & 1)   | (x & 2)) | ((y & 2)   << 1);
	return (val << 2) | (x&1)    | ((y&1)   << 1);

	// 29 shifts, 29 or's, 30 and's
	/*unsigned int val = (x&1) | ((y&1) << 1) | ((z&1) << 2); // non-optimized version
	val |= ((x&2) << 2) | ((y&2) << 3) | ((z&2) << 4);
	val |= ((x&4) << 4) | ((y&4) << 5) | ((z&4) << 6);
	val |= ((x&8) << 6) | ((y&8) << 7) | ((z&8) << 8);
	val |= ((x&16) << 8) | ((y&16) << 9) | ((z&16) << 10);
	val |= ((x&32) << 10) | ((y&32) << 11) | ((z&32) << 12);
	val |= ((x&64) << 12) | ((y&64) << 13) | ((z&64) << 14);
	val |= ((x&128) << 14) | ((y&128) << 15) | ((z&128) << 16);
	val |= ((x&256) << 16) | ((y&256) << 17) | ((z&256) << 18);
	val |= ((x&512) << 18) | ((y&512) << 19) | ((z&512) << 20);
	return val;*/
}
static inline unsigned int unInterleaveBitCompenentX(const unsigned int val)
{
	return (val&0x1)             | ((val&0x8)      >> 2)  | ((val&0x40)      >> 4) | ((val&0x200)     >> 6)  | ((val&0x1000)     >> 8)  |
	       ((val&0x8000)  >> 10) | ((val&0x40000)  >> 12) | ((val&0x200000) >> 14) | ((val&0x1000000) >> 16) | ((val&0x8000000)  >> 18);
}
static inline unsigned int unInterleaveBitCompenentY(const unsigned int val)
{
	return ((val&0x2)     >> 1)  | ((val&0x10)     >> 3)  | ((val&0x80)     >> 5)  | ((val&0x400)     >> 7)  | ((val&0x2000)     >> 9)  |
	       ((val&0x10000) >> 11) | ((val&0x80000)  >> 13) | ((val&0x400000) >> 15) | ((val&0x2000000) >> 17) | ((val&0x10000000) >> 19);
}
static inline unsigned int unInterleaveBitCompenentZ(const unsigned int val)
{
	return ((val&0x4)     >> 2)  | ((val&0x20)     >> 4)  | ((val&0x100)    >> 6)  | ((val&0x800)     >> 8)  | ((val&0x4000)     >> 10) |
	       ((val&0x20000) >> 12) | ((val&0x100000) >> 14) | ((val&0x800000) >> 16) | ((val&0x4000000) >> 18) | ((val&0x20000000) >> 20);
}
static inline void unInterleaveBits(const unsigned int val, unsigned int &x, unsigned int &y, unsigned int &z)
{
	x = unInterleaveBitCompenentX(val);
	y = unInterleaveBitCompenentY(val);
	z = unInterleaveBitCompenentZ(val);
}


static inline double randomNumber(const double lowVal, const double highVal)
{
	return (((double)rand()) / RAND_MAX*(highVal - lowVal) + lowVal);
}
static inline int randomInteger(const int lowVal, const int highVal)
{
	return rand() % (highVal - lowVal + 1) + lowVal;
}


// apparently Visual Studio does not do well with small negative numbers and powers, so use this instead of "cbrt" function
static inline double cubeRoot(const double val)
{
	return (val < 0.0 ? -pow(-val, 1.0 / 3.0) : pow(val, 1.0 / 3.0));
}

// finds real roots of x^3 + A*x^2 + B*x + C
// the values stored in the roots are sorted in descending order
// returns root type
//   CUBIC_ROOT_SINGLE: root0 is real, root1 == root2 == NULL are imaginary
//   CUBIC_ROOT_TRIPLE: root0 == root1 == root2 are real
//   CUBIC_ROOT_SINGLE_DOUBLE: root0, root1 == root2 or root0 == root1, root2 - all are real
//   CUBIC_ROOT_THREE: root0, root1, root2 are all real
// taken from http://teem.sourceforge.net/
#ifndef CUBIC_ROOT_TYPE
#define CUBIC_ROOT_TYPE
enum CubicRootType { CUBIC_ROOT_SINGLE = 0, CUBIC_ROOT_TRIPLE, CUBIC_ROOT_SINGLE_DOUBLE, CUBIC_ROOT_THREE };
#endif

static inline int getCubicRoots(double &root0, double &root1, double &root2, const double A, const double B, const double C,
                                const bool refineWithNewtonRaphson = true, const double eps = 1.0E-11)
{
	double sub = A / 3.0;
	double AA = A * A;
	double Q = (AA / 3.0 - B) / 3.0;
	double R = (-2.0 * A * AA / 27.0 + A * B / 3.0 - C) * 0.5;
	double QQQ = Q * Q * Q;
	double D = R * R - QQQ;

	if (D < -eps) {
		//three distinct roots- this is the most common case, it has been tested the most, its code should go first
		double theta = acos(R / sqrt(QQQ)) / 3.0;
		double t = 2.0 * sqrt(Q);
		// these are sorted, because the C definition of acos says that it returns values in [0, pi]
		root0 = t * cos(theta) - sub;
		root1 = t * cos(theta - M_2_PI / 3.0) - sub;
		root2 = t * cos(theta + M_2_PI / 3.0) - sub;
		return CUBIC_ROOT_THREE;
	}
	else if (D > eps) {
		// one real solution, except maybe also a "rescued" double root
		double sqrt_D = sqrt(D);
		double u = cubeRoot(sqrt_D + R);
		double v = -cubeRoot(sqrt_D - R);
		double x = u + v - sub;

		if (!refineWithNewtonRaphson) {
			root0 = x;
			return CUBIC_ROOT_SINGLE;
		}

		// else refine x, the known root, with newton-raphson to get the most accurate possible calculation for nr, the possible new root
		double der;
		x -= (der = (3.0 * x + 2.0 * A) * x + B, ((x / der + A / der) * x + B / der) * x + C / der);
		x -= (der = (3.0 * x + 2.0 * A) * x + B, ((x / der + A / der) * x + B / der) * x + C / der);
		x -= (der = (3.0 * x + 2.0 * A) * x + B, ((x / der + A / der) * x + B / der) * x + C / der);
		x -= (der = (3.0 * x + 2.0 * A) * x + B, ((x / der + A / der) * x + B / der) * x + C / der);
		x -= (der = (3.0 * x + 2.0 * A) * x + B, ((x / der + A / der) * x + B / der) * x + C / der);
		x -= (der = (3.0 * x + 2.0 * A) * x + B, ((x / der + A / der) * x + B / der) * x + C / der);

		double nr = -(A + x) * 0.5;
		double fnr = ((nr + A) * nr + B) * nr + C;  // the polynomial evaluated at nr

		if (fnr < -eps || fnr > eps) {
			root0 = x;
			return CUBIC_ROOT_SINGLE;
		}
		else {
			if (x > nr) {
				root0 = x;
				root1 = nr;
				root2 = nr;
			}
			else		{
				root0 = nr;
				root1 = nr;
				root2 = x;
			}
			return CUBIC_ROOT_SINGLE_DOUBLE;
		}
	}
	else { // else D is in the interval [-epsilon, +epsilon]
		if (R < -eps || eps < R) { // one double root and one single root
			double u = cubeRoot(R);
			if (u > 0) {
				root0 = 2.0 * u - sub;
				root1 = -u - sub;
				root2 = -u - sub;
			}
			else	   {
				root0 = -u - sub;
				root1 = -u - sub;
				root2 = 2.0 * u - sub;
			}
			return CUBIC_ROOT_SINGLE_DOUBLE;
		}
		else { // one triple root
			root0 = root1 = root2 = -sub;
			return CUBIC_ROOT_TRIPLE;
		}
	}
}
}

#endif


