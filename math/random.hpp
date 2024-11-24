#pragma once

#include <random>
#include <math/types.hpp>

namespace spurt {
    
    template<typename T, typename Enable>
    struct distribution_traits {};
    
    template<typename T>
    struct distribution_traits<T, typename = std::enable_if<std::is_integral<T>::value>::type>
    {
        typedef T value_type;
        typedef std::uniform_int_distribution<T> uniform_distribution;
        static constexpr value_type default_min = 0;
        static constexpr value_type default_max = 100;
    };
    
    
    template<typename T>
    struct distribution_traits<T, typename = std::enable_if<std::is_floating_point<T>::value>::type>
    {
        typedef T value_type;
        typedef std::uniform_real_distribution<T> uniform_distribution;
        static constexpr value_type default_min = 0.;
        static constexpr value_type default_max = 1.;
    };

    
    template<typename T, typename Dist = typename distribution_traits<T>::uniform_distribution>
    class random_generator
    {
    public:
        typedef distribution_traits<T> traits_type;
        typedef T value_type;
        typedef Dist distribution_type;
        
        random_generator(long seed, value_type min = traits_type::default_min, 
                         value_type max = type_traits::default_max) 
            : m_min(min), m_max(max), m_seed(seed), m_random_dev() {
            initialize();
        }
        
        random_generator() : m_min(traits_type::default_min), 
            m_max(traits_type::default_max), m_random_dev() {
            m_seed = m_random_dev();
            initialize();
        }
        
        void initialize() {
            m_random_eng = std::default_random_engine(m_seed);
            m_uniform_dist = std::uniform_int_distribution<value_type>(m_min, m_max);
        }
        
        value_type operator() const {
            return m_uniform_dist(m_random_eng);
        }
        
        
    private:
        std::default_random_engine m_random_eng;
        std::random_device m_random_dev;
        distribution_type m_uniform_dist;
        long m_seed;
        value_type m_min, m_max;
    };
}



























}