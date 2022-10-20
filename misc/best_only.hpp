#ifndef __XAVIER_BEST_ONLY_HPP__
#define __XAVIER_BEST_ONLY_HPP__

namespace spurt {

template<typename Type, typename Value, typename Less=std::less<Value> >
struct best_only {
    best_only(const Less& less = Less()) : __less(less), __empy(true), __t(), __v() {}
    best_only(const best_only& bo) : __less(bo.__less), __empty(bo.__empty),
        __t(bo.__t), __v(bo.__v) {}
        
    void add(const Type& t, const Value& v) {
        if (__empty || __less(t, __t)) {
            __t = t;
            __v = v;
            __empty = false;
        }
    }
    
    std::pair<Type, Value> best() const {
        return std::pair<Type, Value>(__t, __v);
    }
    
    Type    __t;
    Value   __v;
    Less    __less;
};
}

#endif
