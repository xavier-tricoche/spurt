#ifndef __XAVIER_MISC_LOG_HELPER_HPP__
#define __XAVIER_MISC_LOG_HELPER_HPP__

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace spurt { namespace log { 
    
class dual_ostream {
public:
    typedef dual_ostream self_t;
    
    dual_ostream(std::ofstream& logfile, 
                 std::ostream& os=std::cout,
                 int nthreads=1,
                 int os_thresh=0, int file_thresh=1, bool verbose=true) 
        : m_file(logfile), m_os(os), m_file_thresh(file_thresh), 
          m_os_thresh(os_thresh), m_verbose(verbose), m_streams(nthreads),
          m_files(nthreads), m_priority(nthreads), m_active(nthreads, false) {
        m_max_prio = std::max(os_thresh, file_thresh);
        set_nthreads(nthreads);
#ifdef _OPENMP
        omp_init_lock(&m_lock);
#endif
    }
          
    ~dual_ostream() {
        for (int i=0; i<m_streams.size(); ++i) {
            if (m_streams[i]!=NULL) delete m_streams[i];
            if (m_files[i]!=NULL) delete m_files[i];
        }
#ifdef _OPENMP
        omp_destroy_lock(&m_lock);
#endif
    }
    
    void set_nthreads(int n) {
        if (n != m_streams.size()) {
            for (int i=0; i<m_streams.size(); ++i) {
                if (m_streams[i]!=NULL) delete m_streams[i];
                if (m_files[i]!=NULL) delete m_files[i];
            }
            m_streams.resize(n);
            m_priority.resize(n);
            m_files.resize(n);
            m_active.resize(n);
        }
        for (int i=0; i<n; ++i) {
            m_streams[i] = new oss_t();
            m_files[i] = new oss_t();
            m_priority[i]=1000;
            m_active[i]=false;
        }
    }
          
    void set_thresholds(int os_thresh, int file_thresh) {
        m_os_thresh = os_thresh;
        m_file_thresh = file_thresh;
        m_max_prio = std::max(os_thresh, file_thresh);
        for (int i=0; i<m_active.size(); ++i) {
            m_active[i]=false;
        }
    }
    
    self_t& operator()(int p=0, const std::string& pre="") {
        const int thread=get_id();
        m_priority[thread]=p;
        if (p > m_max_prio) {
            m_active[thread] = false;
            return *this;
        }
        else m_active[thread] = true;
        
        m_streams[thread]->clear();
        m_files[thread]->clear();
        if (p <= m_os_thresh) {
            *(m_streams[thread]) << pre << "[thread " << thread << "]: ";
        }
        if (p <= m_file_thresh) {
            *(m_files[thread]) << pre << "[thread " << thread << "]: ";
        }
        return *this;
    }
    
    self_t& operator<<(std::ios_base& (*func)(std::ios_base&) ) {
        const int thread = get_id();
        
        if (!m_active[thread]) return *this;
        if (m_priority[thread] <= m_os_thresh) {
            func(*m_streams[thread]);
        }
        if (m_priority[thread] <= m_file_thresh) {
            func(*m_files[thread]);
        }
        return *this;
    }
    
    template<typename CharT=char, typename Traits=std::char_traits<CharT> > 
    self_t& operator<<(std::basic_ios<CharT,Traits>& (*func)(std::basic_ios<CharT,Traits>&) ) {
        const int thread = get_id();
        
        if (!m_active[thread]) return *this;
        if (m_priority[thread] <= m_os_thresh) {
            func(*m_streams[thread]);
        }
        if (m_priority[thread] <= m_file_thresh) {
            func(*m_files[thread]);
        }
        return *this;
    }
    
    template<typename CharT=char, typename Traits=std::char_traits<CharT> > 
    self_t& operator<<(std::basic_ostream<CharT,Traits>& (*func)(std::basic_ostream<CharT,Traits>&) ) {
        const int thread = get_id();
        
        if (!m_active[thread]) return *this;
        if (m_priority[thread] <= m_os_thresh) {
            func(*m_streams[thread]);
#ifdef _OPENMP
            omp_set_lock(&m_lock);
#endif
            m_os << m_streams[thread]->str() << std::flush;
#ifdef _OPENMP
            omp_unset_lock(&m_lock);
#endif
        }
        if (m_priority[thread] <= m_file_thresh) {
#ifdef _OPENMP
	    omp_set_lock(&m_lock);
#endif
            func(*m_files[thread]);
            m_file << m_files[thread]->str() << std::flush;
#ifdef _OPENMP
            omp_unset_lock(&m_lock);
#endif
        }
        m_active[thread] = false;
        m_streams[thread]->str("");
        m_streams[thread]->clear();
        m_files[thread]->str("");
        m_files[thread]->clear(); 
        return *this;
    }
    
    template< typename T >
    friend self_t& operator<<(self_t&, const T&);
    
private:
    typedef std::ostringstream oss_t;
    
#ifdef _OPENMP
    omp_lock_t m_lock;
#endif
    
    int m_file_thresh;
    int m_os_thresh;
    int m_p;
    int m_max_prio;
    bool m_verbose;
    
    std::ostream& m_os;
    std::ofstream& m_file;
    
    // per-thread storage
    std::vector< oss_t* > m_streams;
    std::vector< oss_t* > m_files;
    std::vector< int > m_priority;
    std::vector< bool > m_active;
    
    void reset() {
        m_p=0;
    }
    
    int get_id() const {
#if _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }
};

template< typename T >
dual_ostream& operator<<(dual_ostream& os, const T& t) {
#if _OPENMP
        const int thread = omp_get_thread_num();
#else
        const int thread = 0;
#endif
    if (os.m_priority[thread] <= os.m_os_thresh) {
        *(os.m_streams[thread]) << t;
    }
    if (os.m_priority[thread] <= os.m_file_thresh) {
        *(os.m_files[thread]) << t;
    }
    return os;
}


} // log
} // xavier



#endif

