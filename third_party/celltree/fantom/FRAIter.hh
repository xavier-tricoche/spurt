template<class RAContainer,class T>
class CIT{

protected:

  RAContainer * that; ///< this is not const due to our derived class IT
  positive index;   

public:

  typedef T value_type;
  typedef T* pointer;
  typedef T& reference;
  typedef const T& const_reference;
  typedef positive size_type;
  typedef long difference_type;
  typedef random_access_iterator_tag iterator_category;
  

  /// constructor takes a const pointer, which is converted to a non-const
  /// pointer if the iterator is not a const iterator
  /// the problem here is, that we might call the non-const version of a
  /// function which might have a different semantic than the const version
  inline CIT(const RAContainer * ref=0,positive ind = 0)
    :that(const_cast<RAContainer*>(ref)),index(ind)   
  {}    

  // access functions

  inline const_reference operator * () const { return (*that)[index]; }

  inline const pointer operator ->() const { return &((*that)[index]); }

  inline const_reference operator [] (positive i) const { return (*that)[index P i]; }


  // iteration functions

  inline CIT & operator ++()  { PP index;return *this;} 
  inline CIT & operator --()  { MM index;return *this;}

  inline CIT operator ++( int ) { CIT oit = *this;  PP index; return oit;  }
  inline CIT operator --( int ) { CIT oit = *this;  MM index; return oit;  }

  inline CIT operator + (positive b)  { return CIT(that,index P b); }
  inline CIT operator - (positive b)  { return CIT(that,index M b); }

  inline CIT operator += (positive b) { index PE b; return *this; }
  inline CIT operator -= (positive b) { index ME b; return *this; }    


  //comparisons

  inline bool operator != (const CIT & b) const { return index != b.index  |  that != b.that; }
  inline bool operator <  (const CIT & b) const { return index+1 LT b.index+1; } // +1 on both sides to ensure that -1 as rend() is handled correctly
  inline bool operator == (const CIT & b) const { return index == b.index  &  that == b.that; }

  //difference
  inline difference_type operator-(const CIT & b) { return P(index-b.index); }

  friend difference_type operator-( const CIT<RAContainer,T> & a, const CIT<RAContainer,T> & b )
  {
      return P( a.index - b.index );
  }
};

template<class RAContainer,class T>
class IT:public CIT<RAContainer,T> {

public:

  typedef random_access_iterator_tag iterator_category;

  inline IT(RAContainer * ref=0,positive ind = 0)
    :CIT<RAContainer,T>(ref,ind)   
  {}    

  inline T & operator * () { return (*(this->that))[this->index]; }
  inline T & operator ->() { return &((*(this->that))[this->index]); }
  inline T & operator [] (positive i) { return (*(this->that))[this->index P i]; }
    
  inline IT & operator ++()  { PP this->index;return *this;} 
  inline IT & operator --()  { MM this->index;return *this;}

  inline IT operator ++( int ) { IT oit = *this;  PP this->index; return oit;  }
  inline IT operator --( int ) { IT oit = *this;  MM this->index; return oit;  }

  inline IT operator + (positive b)  { return IT(this->that,this->index P b); }
  inline IT operator - (positive b)  { return IT(this->that,this->index M b); }

  inline IT operator += (positive b) { this->index PE b; return *this; }
  inline IT operator -= (positive b) { this->index ME b; return *this; }    


};


