#ifndef F_ANY_ARRAY_HH
#define F_ANY_ARRAY_HH

#include<vector>
#include "stdAliases.hh"
#include<iterator>
#include "FException.hh"
using namespace std;
#include <iostream> // debug only


namespace FIter{

#define P +
#define M -
#define PE +=
#define ME -=
#define PP ++
#define MM --
#define LT < 
#define IT iterator 
#define CIT const_iterator 

#include "FRAIter.hh"

#undef P
#undef M
#undef PE
#undef ME
#undef PP
#undef MM
#undef LT
#undef IT
#undef CIT

#define P -
#define M +
#define PE -=
#define ME +=
#define PP --
#define MM ++
#define LT > 
#define IT reverse_iterator 
#define CIT const_reverse_iterator 

#include "FRAIter.hh"

#undef P
#undef M
#undef PE
#undef ME
#undef PP
#undef MM
#undef LT
#undef IT
#undef CIT
  


}



template<class T> 
class FanyArray
{

 public:
  typedef T value_type;
  typedef T* pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef positive size_type;
  typedef long difference_type;

  typedef FIter::const_iterator<FanyArray<T>,T> const_iterator;

  typedef FIter::const_reverse_iterator<FanyArray<T>,T> const_reverse_iterator;

  typedef positive distance_type;

  //virtual T & operator [] ( positive ) = 0;
  //virtual const T & operator [] ( positive ) const = 0;

  virtual const T operator [] ( positive ) const = 0;

  virtual const T & operator () ( positive i ) const 
  {
    static T x;
    x=(*this)[i];
    return x;
  };
  
  //copies n values of type t to memory area pointed to by buf
  virtual void get_range ( positive firstIndex, positive nIndices, T*buf ) const
  { for(positive i=0;i<nIndices;i++) buf[i] = (*this)[firstIndex+i]; } 

  virtual positive size() const =0;
  virtual positive max_size() const =0;

  inline T front() const {return (*this)[0]; }
  inline T back() const {return (*this)[size()-1]; }
          		       
  inline const_iterator begin() const { return const_iterator(this,0);    }
  inline const_iterator end() const { return const_iterator(this,size()); }   

  inline const_reverse_iterator rbegin() const { return const_reverse_iterator(this,0);    }
  inline const_reverse_iterator rend() const { return const_reverse_iterator(this,size()); }   

  virtual positive memSize() const { std::cout << "FanyArray::memSize(): oops, this should be pure virtual:\n"; return 0; } // dummy, should be pure virtual
  virtual ~ FanyArray() {};

};

template<class T>
class FAssignableAnyArray : public FanyArray<T>
{
  public:
  typedef FIter::iterator<FAssignableAnyArray<T>,T> iterator;
  typedef FIter::reverse_iterator<FAssignableAnyArray<T>,T> reverse_iterator;

  inline iterator begin() { return iterator(this,0);      }
  inline iterator end() { return iterator(this,this->size()); }   

  inline reverse_iterator rbegin() { return reverse_iterator(this,0);      }
  inline reverse_iterator rend() { return reverse_iterator(this,this->size()); }   

  inline T& front() {return (*this)[0]; }
  inline T& back() {return (*this)[this->size()-1]; }

  virtual T & operator [] ( positive ) = 0;
  virtual positive memSize() const { return FanyArray<T>::memSize(); } // dummy
  virtual ~FAssignableAnyArray(){}
};


template<class T>
class FanyVector:public FAssignableAnyArray<T>
{
  vector<T> vec;
  public:
  FanyVector(){}
  FanyVector(vector<T> & v) {vec.swap(v);}
  FanyVector(positive n):vec(n){ }

  T & operator[] ( positive i )
  {
#ifndef NODEBUG
    if(i>=vec.size()){
      THROW_DEFAULT_EXCEPTION(FIndexOutOfBoundsException);
    }
#endif
    return vec[i];
  }

  //const T & operator [] ( positive i ) const
  const T operator [] ( positive i ) const
  {

#ifndef NODEBUG
    if(i>=vec.size()){
      THROW_DEFAULT_EXCEPTION(FIndexOutOfBoundsException);
    }
#endif
    return vec[i];
  }

  //const T & operator [] ( positive i ) const
  const T & operator () ( positive i ) const
  {

#ifndef NODEBUG
    if(i>=vec.size()){
      THROW_DEFAULT_EXCEPTION(FIndexOutOfBoundsException);
    }
#endif
    return vec[i];
  }

  positive size() const
  {return vec.size();}

  virtual void get_range(positive i,positive n,T*buf) const
  {
    memcpy(buf,&vec[i],sizeof(T)*n);
  }

  positive max_size() const
  {return vec.size();}

  positive memSize() const
  { return vec.size() * sizeof(T); }

  virtual ~FanyVector(){}
};

#endif //F_ANY_ARRAY_HH
