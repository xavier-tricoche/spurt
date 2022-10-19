#ifndef F_ANY_MMAPPED_ARRAY_HH
#define F_ANY_MMAPPED_ARRAY_HH

#include "FanyArray.hh"
#include "FmmappedArray.hh"
#include <sys/stat.h>

template<class T>
class FanyMMappedArray:public FAssignableAnyArray<T>
{

public:

  mutable shared_ptr<FmmappedArray> arr;

  FanyMMappedArray()
  {
  }

  FanyMMappedArray(const string& name,bool readonly=true)
  {
    struct stat s;
    stat(name.c_str(),&s);              
    arr.reset(new FmmappedArray(vector<string>(1,name),
				vector<positive>(1,s.st_size/sizeof(T)),
				sizeof(T),128<<10,10,readonly));
  }

  FanyMMappedArray(vector<string> names,vector<positive> sizes
		   =vector<positive>(),bool readonly=false)
  {
    
    if(sizes.empty()){
      
      sizes.resize(names.size());

      for(unsigned i=0;i<names.size();i++){
	struct stat s;
	stat(names[i].c_str(),&s);      
	sizes[i] = s.st_size/sizeof(T);
	
      }		
      
    }
    

    //erase 0-sized files from lists

    size_t writei = 0;
    for(size_t i=0;i<sizes.size();++i)
      if( sizes[i]!=0 ){ 
	if(writei < i )
	  {
	    sizes[writei]=sizes[i];
	    names[writei]=names[i];
	  }
	++writei;
      }
    names.resize(writei);
    sizes.resize(writei);
	

      

    if(readonly)
      arr.reset(new FmmappedArray(names,sizes,sizeof(T),512<<10,10,readonly));
    else
      arr.reset(new FmmappedArray(names,sizes,sizeof(T)));

  }

  FanyMMappedArray(shared_ptr<FmmappedArray> a)
    :arr(a)
  {
    assert(a->getBlockSize()==sizeof(T));
  }

  const T operator [] ( positive i ) const
  { 
    //!! attention : if this operator is invoked a 2nd time, it is 
    //  __not guaranteed__ that the reference stays valid !! 
    positive n;
    return *((const T *)arr->getBlock(n,i));
  }

  const T & operator() ( positive i )
  { 
    //!! attention : if this operator is invoked a 2nd time, it is 
    //  __not guaranteed__ that the reference stays valid !! 
    positive n;
    return *((T*)arr->getBlock(n,i));
  }

  T & operator[] ( positive i )
  { 
    positive n;
    return *((T*)arr->getBlock(n,i));
  }

  virtual void get_range ( positive firstIndex, positive nIndices, T*buf ) const
  {
    positive n;
    void*p;
    for(;;){
      p = (T*)arr->getBlock(n,firstIndex);
      if(n>=nIndices){
	memcpy(buf,p,nIndices*sizeof(T));
	return;
      }
      memcpy(buf,p,n*sizeof(T));
      firstIndex+=n;
      buf+=n;
      nIndices-=n;
    }    
  }

  positive size() const
  { return arr->getNbBlocks(); }

  positive max_size() const
  { return arr->getNbBlocks(); }

  ~FanyMMappedArray(){}

};

#endif //F_ANY_MMAPPED_ARRAY_HH
