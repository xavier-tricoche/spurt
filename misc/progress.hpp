#ifndef __PROGRESS_HPP__
#define __PROGRESS_HPP__

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <chrono>

namespace spurt {
    
class timer {
    typedef std::chrono::high_resolution_clock hr_clock_t;
    typedef hr_clock_t::time_point time_point_t;
    
    std::clock_t m_cpu_begin, m_cpu_end; 
    time_point_t m_wall_begin, m_wall_end;
    
    template<typename _Time>
    double __elapsed(const _Time& begin, const _Time& end) {
        return 1000.*(end-begin)/CLOCKS_PER_SEC;
    }
    
    bool m_stopped;
        
public:
    timer() : m_stopped(false) {}
    
    time_point_t wall_now() const {
        return hr_clock_t::now();
    }
    
    std::clock_t cpu_now() const { 
        return std::clock();
    }
    
    void start() {
        m_cpu_begin = cpu_now();
        m_wall_begin = wall_now();
        m_stopped = false;
    }
    
    void stop() {
        m_cpu_end = cpu_now();
        m_wall_end = wall_now();
        m_stopped = true;
    }
    
    double cpu_elapsed() const {
        auto cur = m_stopped ? m_cpu_end : cpu_now();
        return elapsed(m_cpu_begin, cur);
    }
    
    double wall_elapsed() const {
        auto cur = m_stopped ? m_wall_end : wall_now();
        return elapsed(m_wall_begin, cur);
    }
    
    double elapsed() const { 
        return wall_elapsed();
    }
};

struct ProgressDisplay {
    typedef std::chrono::high_resolution_clock hr_clock_t;
    typedef hr_clock_t::time_point time_point_t;
    
    bool active;
    float progress;
    int bar_width;
    int precision;
    size_t size;
    std::clock_t cpu_time_begin, cpu_time_end; 
    time_point_t wall_time_begin, wall_time_end;

    std::ostream& m_os;
    std::string m_str;
    int m_pos;
    size_t m_at, m_delta_at;
      
    ProgressDisplay(bool activate=true, std::ostream& _os=std::cout) 
        : active(activate), progress(0), bar_width(60), precision(3), 
          m_os(_os), m_str("") {}

    void start(size_t _size, const std::string& str="", size_t nb_updates=1000) {
        m_str=str;
        size=_size;
        progress=0;
        m_at = 0;
        m_delta_at = std::max(static_cast<size_t>(1), _size/nb_updates);
        cpu_time_begin = std::clock();
        wall_time_begin = std::chrono::high_resolution_clock::now();
        if (active) display(); 
    }

    void end() { 
        cpu_time_end = std::clock();
        wall_time_end = std::chrono::high_resolution_clock::now();
        if (active) {
            display();
            std::cout 
                << "\ncompute time: cpu: " 
                << cpu_time() 
                << " ms. (" << 1000.*static_cast<double>(size)/cpu_time() 
                << " Hz) | wall: " 
                << wall_time()
                << " ms." << std::endl;
        }
    }

    void update(size_t at) {
        if (at - m_at >= m_delta_at) {
            progress=(float)(at+1)/(float)size;
            m_at = at;
            if (active) display(); 
        }
    }
    
    double cpu_time() const {
        return 1000.*(cpu_time_end-cpu_time_begin)/CLOCKS_PER_SEC;
    }
    
    double instant_cpu_time() const {
        return 1000.*(std::clock()-cpu_time_begin)/CLOCKS_PER_SEC;
    }
    
    double wall_time() const {
        return std::chrono::duration<double, std::milli>(wall_time_end-wall_time_begin).count();
    }
    
    double instant_wall_time() const {
        return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-wall_time_begin).count();
    }
    
    void display() {
        if (!m_str.empty()) m_os << m_str << ": ";
        m_os << "[";
        int pos=bar_width * progress;
        m_os << std::string(pos, '=') << '>' 
           << std::string(std::max(bar_width-pos-1, 0), ' ')
           << "] " 
           << std::setprecision(precision)
           << std::fixed
           << progress*100.0
           << " %\r"
           << std::flush;
    }
};

} // namespace spurt

#endif
