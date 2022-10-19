#ifndef __XAVIER_BEST_ONLY_HPP__
#define __XAVIER_BEST_ONLY_HPP__

namespace spurt {

template<typename Type, typename Value, typename Less=std::less<Value> >
struct best_only {
    best_only(const Less& less = Less()) : m_less(less), m_empy(true), m_t(), m_v() {}
    best_only(const best_only& bo) : m_less(bo.m_less), m_empty(bo.m_empty),
        m_t(bo.m_t), m_v(bo.m_v) {}
        
    void add(const Type& t, const Value& v) {
        if (m_empty || m_less(t, m_t)) {
            m_t = t;
            m_v = v;
            m_empty = false;
        }
    }
    
    std::pair<Type, Value> best() const {
        return std::pair<Type, Value>(m_t, m_v);
    }
    
    Type    m_t;
    Value   m_v;
    Less    m_less;
};
}

#endif
