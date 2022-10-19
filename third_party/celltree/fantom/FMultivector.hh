//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FMultivector.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:27 $
// Author:    $Author: garth $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#ifndef __FMultivector_hh
#define __FMultivector_hh

#include "stdAliases.hh"
#include "FArray.hh"
#include "FTensor.hh"

#include <vector>

class FRefTensor;

//===========================================================================

#include "FVector.hh"
// class FMatrix;

/** 
 *The FMultivector class provides a representation of the mathematical
 *multivector. Its order is 1, its dimension 4 or 8;
 */
class FMultivector : public FTensor
{
    public:

        /** 
         *\par Description:
         *Constructor: provides an empty multivector.
         */
        FMultivector();

        /** 
         *\par Description:
         *Constructor: provides a multivector of dimension \b dimension and order
         *\b order.
         *\pre
         *order ==1, dimension = {4,8}
         *\post
         *none
         *\exception
         *FInvalidDimensionException
         *\param
         *dimension: multivector dimension
         *\param
         *order: multivector order
         *\param
         *clear: boolean to indicate that the constructor must return a zero multivector
         */
        FMultivector(unsigned char dimension, unsigned char order, bool clear=false);

        /** 
         *\par Description:
         *Constructor: provides a multivector of dimension \b dimension and order 
         *\b order and with double coordinates \b comp.
         *\pre
         *order ==1, dimension = {4,8}\\
         *The size of comp must fit the one defined by order and dimension
         *\post
         *none
         *\exception
         *FInvalidDimensionException
         *\param
         *dimension: value of the dimension.
         *\param
         *order: value of the order.
         *\param
         *comp: array of double coordinates.
         */
        FMultivector(unsigned char dimension, unsigned char order, 
                const std::vector<double>& a);

        /** 
         *\par Description:
         *Copy constructor: provides a multivector, copy of its argument. 
         *\pre
         *none
         *\post
         *none
         *\exception
         *none
         *\param
         *T: multivector to copy.
         */
        FMultivector(const FMultivector& T);

        /**
         * Copy constructor. Generate multivector from tensor.
         * \param m
         * Tensor to copy.
         */
        explicit FMultivector(const FTensor& m);

        /**
         * Copy constructor. Generate multivector from vector.
         * \param v
         * Vector to copy.
         */
        explicit FMultivector(const FVector& v);

        /** 
         *\par Description:
         *Destructor.
         *\pre
         *none
         *\post
         *none
         *\exception
         *none
         */
        ~FMultivector();

        /** 
         *\par Description:
         *Gets the class name.
         *\pre
         *none
         *\post
         *none
         *\exception
         *none
         */
        //  virtual const FString& getClassName() const; 


        /** 
         * sets dimension of the multivector (and invalidates the actual value !)
         * \exception
         * FInvalidDimensionException, FException...
         */
        void setDimension(unsigned char );

        /** 
         * Sets order of the multivector (and invalidates the actual value !)
         * \exception
         * FInvalidDimensionException, FException...
         */
        void setOrder(unsigned char );

        /** resizer(version of resize ( , ) for compatibility.
         * \param dim
         * New dimension.
         * \param ord
         * New order.
         */
        void resizeMultivector (unsigned char dim, unsigned char ord);


        /** 
         *\par Description:
         *Assignment of a multivector.
         *\pre
         *none
         *\post
         *Every components have been copied.
         *\exception
         *none
         *\param
         *T: tensor to copy.
         */
        FMultivector& operator=(const FMultivector& T);

        /** 
         *\par Description:
         *Set all components of a multivector to a given double value
         *\pre
         *multivector's dimension and order have been set 
         *\post
         *Every components have been set
         *\exception
         *FEmptyObjectException
         *\param
         *val: double value.
         */
        FMultivector& operator=(double val);

        /** 
         *\par Description:
         * Dummy
         */
        const double& operator()(void) const; 

        /** 
         *\par Description:
         *Gets a scalar component of a multivector.
         *\pre
         *the multivector is of order 1.\\
         *i < dimension
         *\post
         *none
         *\exception
         *FInvalidDimensionException
         *\exception
         *FInvalidIteratorException
         *\exception
         *FEmptyObjectException
         *\param
         *i: index of the component to get.
         */
        const double& operator()(unsigned char i) const; 
        double& operator()(unsigned char i); 

        /** 
         *\par Description:
         * Dummy
         */
        const double& operator()(unsigned char i, unsigned char j) const; 

        /** 
         *\par Description:
         * Dummy
         */
        const double& operator()(unsigned char i, unsigned char j, unsigned char k) const; 


        /** 
         *\par Description:
         * Dummy
         */
        const FRefTensor operator [] (unsigned char i) const;

        /** 
         *\par Description:
         * Dummy
         */
        FRefTensor operator [] (unsigned char i);



        /** 
         *\par Description:
         *prints the contents of the position to os
         *\pre
         *none
         *\post
         *none
         *\exception
         *none
         *\param multivector to print.
         */
        friend std::ostream& operator<< (std::ostream &os, const FMultivector &tensor);

        /** 
         *\par Description:
         * Dummy
         */
        void setValue(double val); 

        /** 
         *\par Description:
         *Sets a component of a multivector.
         *\pre
         *the multivector is of order 1
         *\post
         *none
         *\exception
         *FInvalidDimensionException
         *\exception
         *FInvalidIteratorException
         *\exception
         *FEmptyObjectException
         *\param
         *i: index of the component to get.
         *\param
         *val: scalar value to set.
         */
        void setValue(unsigned char i, double val); 

        /** 
         * Dummy
         */
        void setValue(unsigned char i, unsigned char j, double val); 

        /** 
         *\par Description:
         * Dummy
         */
        void setValue(unsigned char i, unsigned char j, unsigned char k, double val); 

        /** 
         *\par Description:
         *Sets all the scalar components of the multivector.
         *\pre
         *\b comp corresponds to the order and dimension of the multivector.
         *\post
         *none
         *\exception
         *FInvalidDimensionException
         *\param
         *comp: double vector to set all the components of the tensor.
         */
        void setValues(const std::vector<double>& comp);

        /**
         *	\par Description:
         *	Clifford Multiplication
         */
        static FMultivector cliffordMult(const FMultivector& A, const FMultivector& B);

        /**
         *	\par Description:
         *	conjugates the multivector
         */	
        void conjugate();

        /**
         *	\par Description:
         *	returns scalar part of the multivector
         */
        double getScalar() const;

        /**
         *	\par Description:
         *	returns vector part of the multivector
         */
        FVector getVector() const;

        /**
         *	\par Description:
         *	returns bivector part of the multivector
         */
        double getBivector2D() const;

        /**
         *	\par Description:
         *	returns bivector part of the multivector
         */	
        FVector getBivector3D() const;

        /**
         *	\par Description:
         *	returns trivector part of the multivector
         */	
        double getTrivector() const; 

        /// Undocumented.
        friend class FRefTensor;

        /// ScalarProduct
        friend double operator*(const FMultivector& A, const FMultivector& B); 

        /** 
         *\par Description:
         * Dummy 
         */
        friend FMultivector deviator(const FMultivector& tensor);

    private:
        //to prevent invocation of resize of superclass
        void resize(unsigned int s, bool keepContent=false);	 
        friend std::ostream & binwrite( std::ostream & os, const FMultivector & xr );
        friend std::istream & binread( std::istream & is, FTensor * & xpr );
};

//#define __FMultivector_hh_defined

#endif // __FMultivector_hh

//===========================================================================
#ifndef __FMultivector_hh_without_icc

#ifndef OUTLINE
#include "FMultivector.icc"
#endif

#endif  // __FMultivector_hh_without_icc
//=========================================================================== 

