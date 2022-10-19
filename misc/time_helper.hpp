#ifndef __TIME_HELPER_HPP__
#define __TIME_HELPER_HPP__

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
    
    double __elapsed(const time_point_t& begin, const time_point_t& end) const {
        return  std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
    }
    
    double __elapsed(const std::clock_t& begin, const std::clock_t& end) const {
        return (end-begin)/CLOCKS_PER_SEC;
    }
    
    time_point_t __wall_now() const {
        return hr_clock_t::now();
    }
    
    std::clock_t __cpu_now() const { 
        return std::clock();
    }
    
    bool m_stopped;
    
public:
    timer() : m_stopped(false) {
        this->start(); // start timer by default
    }
        
    static inline std::string now() {
        std::time_t t = std::time(nullptr);
        return std::asctime(std::localtime(&t));
    }
    
    void start() {
        m_cpu_begin = __cpu_now();
        m_wall_begin = __wall_now();
        m_stopped = false;
    }
    
    void restart() {
        start();
    }
    
    void stop() {
        m_cpu_end = __cpu_now();
        m_wall_end = __wall_now();
        m_stopped = true;
    }
    
    double cpu_elapsed() const {
        auto cur = m_stopped ? m_cpu_end : __cpu_now();
        return __elapsed(m_cpu_begin, cur);
    }
    
    double wall_elapsed() const {
        auto cur = m_stopped ? m_wall_end : __wall_now();
        return __elapsed(m_wall_begin, cur);
    }
    
    double elapsed() const { 
        return wall_elapsed();
    }
};

struct progress_display : public timer {
    typedef timer base_type;
    
    bool   m_display_on;
    float  m_progress;
    int    m_bar_width;
    int    m_precision;
    size_t m_size;

    std::ostream& m_os;
    std::string   m_str;
    int           m_pos;
    size_t        m_at, m_delta_at;
      
    progress_display(bool display=true, std::ostream& os=std::cout) 
        : base_type(), m_display_on(display), m_progress(0), 
          m_bar_width(60), m_precision(3), m_os(os), m_str("") {}

    void start(size_t size, const std::string& str="", size_t nb_updates=1000) {
        m_str=str;
        m_size=size;
        m_progress=0;
        m_at = 0;
        m_delta_at = std::max(static_cast<size_t>(1), size/nb_updates);
        base_type::start();
        if (m_display_on) display(); 
    }

    void stop() { 
        base_type::stop();
        if (m_display_on) {
            display();
            double cpu_time = base_type::cpu_elapsed();
            double wall_time = base_type::wall_elapsed();
            m_os 
                << "\ncompute time: cpu: " 
                << cpu_time 
                << " ms. (" << 1000.*static_cast<double>(m_size)/cpu_time 
                << " Hz) | wall: " 
                << wall_time
                << " ms." << std::endl;
        }
    }

    void update(size_t at) {
        if (at - m_at >= m_delta_at) {
            m_progress=(float)(at+1)/(float)m_size;
            m_at = at;
            if (m_display_on) display(); 
        }
    }
    
    double cpu_time() const {
        return base_type::cpu_elapsed();
    }
    
    double wall_time() const {
        return base_type::wall_elapsed();
    }
    
    void display() {
        if (!m_str.empty()) m_os << m_str << ": ";
        m_os << "[";
        int pos=m_bar_width * m_progress;
        m_os << std::string(pos, '=') << '>' 
           << std::string(std::max(m_bar_width-pos-1, 0), ' ')
           << "] " 
           << std::setprecision(m_precision)
           << std::fixed
           << m_progress*100.0
           << " %\r"
           << std::flush;
    }
};

} // namespace spurt

#endif
