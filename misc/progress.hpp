#ifndef __PROGRESS_HPP__
#define __PROGRESS_HPP__

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>

#ifdef __OLD_C11_COMPILER__
// missing put_time in standard library
namespace std {
template<typename CharT>
struct _put_time
{
  const std::tm* time;
  const char *fmt;
};
template<typename CharT>
inline _put_time<CharT> put_time(const std::tm* c_time, const CharT* fmt) {
    return { c_time, fmt };
}
template<typename CharT>
std::ostream& operator<<(std::ostream& os, _put_time<CharT> f) {
    char buffer[24];
    if (strftime(buffer, sizeof(buffer), f.fmt, f.time) > 0) {
        os << buffer;
    }
    return os;
}
}
#endif

namespace xavier {

static std::string human_readable_duration(double _t) {
    using namespace std::chrono;
    using dms = duration< double, std::ratio<1, 1000> >;

    if (_t < 0.001) return std::string("0 ms.");
    dms t(_t);
    std::ostringstream os;
    size_t h = duration_cast<hours>(t).count();
    if (h>0) {
        os << h << " h. ";
        t -= duration_cast<milliseconds>(hours(h));
    }
    size_t m = duration_cast<minutes>(t).count();
    if (m>0) {
        os << m << " m. ";
        t -= duration_cast<milliseconds>(minutes(m));
    }
    size_t s = duration_cast<seconds>(t).count();
    if (s>0) {
        os << s << " s. ";
        t -= duration_cast<milliseconds>(seconds(s));
    }
    size_t ms = duration_cast<milliseconds>(t).count();
    if (ms>0) {
        os << ms << " ms. ";
    }
    std::string str=os.str();
    return str.substr(0, str.size()-1);
}

static std::string current_date_and_time(const std::time_t t = std::time(nullptr)) {
    std::ostringstream os;
    os << "UTC    : " << std::put_time(std::gmtime(&t), "%c %Z") << '\n';
    os << "local  : " << std::put_time(std::localtime(&t), "%c %Z") << '\n';
    return os.str();
}

struct ProgressDisplay {
    typedef std::chrono::high_resolution_clock hr_clock_t;
    typedef hr_clock_t::time_point time_point_t;

    bool m_active, m_fraction, m_stopped;
    float m_progress;
    int m_bar_width;
    int m_precision;
    size_t m_size;
    std::clock_t m_cpu_time_begin, m_cpu_time_end;
    time_point_t m_wall_time_begin, m_wall_time_end;
    std::tm *m_begin_date, *m_end_date;

    std::ostream& m_os;
    std::string m_str;
    int m_pos;
    size_t m_at, m_delta_at;

    ProgressDisplay(bool activate=true, std::ostream& _os=std::cout)
        : m_active(activate), m_progress(0), m_bar_width(60), m_precision(3),
          m_os(_os), m_str(""), m_fraction(false), m_stopped(false) {
              m_cpu_time_begin = std::clock();
              m_wall_time_begin = std::chrono::high_resolution_clock::now();
              std::time_t now = std::time(nullptr);
              m_begin_date = std::gmtime(&now);
          }

    void begin(size_t size=0, const std::string& str="", size_t nb_updates=1000) {
        m_str=str;
        m_size=size;
        m_progress=0;
        m_at = 0;
        m_stopped = false;
        m_delta_at = std::max(static_cast<size_t>(1), m_size/nb_updates);
        m_cpu_time_begin = std::clock();
        m_wall_time_begin = std::chrono::high_resolution_clock::now();
        std::time_t now = std::time(nullptr);
        m_begin_date = std::gmtime(&now);
        if (m_active) display();
    }
    void start(size_t size=0, const std::string& str="", size_t nb_updates=1000) {
        this->begin(size, str, nb_updates);
    }

    void end() {
        m_cpu_time_end = std::clock();
        m_wall_time_end = std::chrono::high_resolution_clock::now();
        std::time_t now = std::time(nullptr);
        m_end_date = std::gmtime(&now);
        m_stopped = true;
        if (m_active) {
            display();
            std::cout
                << "\ncompute time: cpu: "
                << human_readable_duration(cpu_time())
                << " (" << 1000.*static_cast<double>(m_size)/cpu_time()
                << " Hz) | wall: "
                << human_readable_duration(wall_time()) << '\n'
                << "Computation started " << std::put_time(m_begin_date, "%c %Z") << ", ended " << std::put_time(m_end_date, "%c %Z") << std::endl;
        }
    }
    void stop() { this->end(); }

    size_t size() const { return m_size; }

    void fraction_on() {
        m_fraction = true;
    }

    void update(size_t at) {
        if (at - m_at >= m_delta_at || at+1 == m_size) {
            m_progress=(float)(at+1)/(float)m_size;
            m_at = at;
            if (m_active) display();
        }
    }

    void set_active(bool active=true) { m_active = active; }

    std::tm* begin_date() const {
        return m_begin_date;
    }

    std::tm* end_date() const {
        if (!m_stopped) {
            std::time_t tmp = std::time(nullptr);
            return std::gmtime(&tmp);
        }
        else return m_end_date;
    }

    double cpu_time() const {
        if (!m_stopped) return instant_cpu_time();
        return 1000.*(m_cpu_time_end-m_cpu_time_begin)/CLOCKS_PER_SEC;
    }

    double instant_cpu_time() const {
        return 1000.*(std::clock()-m_cpu_time_begin)/CLOCKS_PER_SEC;
    }

    double wall_time() const {
        if (!m_stopped) return instant_wall_time();
        return std::chrono::duration<double, std::milli>(m_wall_time_end-m_wall_time_begin).count();
    }

    double instant_wall_time() const {
        return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now()-m_wall_time_begin).count();
    }

    void display() {
        std::string optional;
        if (m_fraction) {
            optional = " (" + std::to_string(m_at) + " / " + std::to_string(m_size) + ")";
        }
        if (!m_str.empty()) m_os << m_str << ": ";
        m_os << "[";
        int pos=m_bar_width * m_progress;
        m_os << std::string(pos, '=') << '>'
           << std::string(std::max(m_bar_width-pos-1, 0), ' ')
           << "] "
           << std::setprecision(m_precision)
           << std::fixed
           << m_progress*100.0
           << " %" << optional << "\r"
           << std::flush;
    }
};

std::ostream& operator<<(std::ostream& os, const ProgressDisplay& pd) {
    os << "wall time: " << human_readable_duration(pd.wall_time())
    << ", cpu time: " << human_readable_duration(pd.cpu_time());
    if (pd.size()>0)
        os << " (" << 1000.*static_cast<double>(pd.size())/pd.cpu_time() << " Hz)";
    os << ", started " << std::put_time(pd.begin_date(), "%c %Z") << ", ended "
    << std::put_time(pd.end_date(), "%c %Z");
    return os;
}

} // namespace xavier

#endif
