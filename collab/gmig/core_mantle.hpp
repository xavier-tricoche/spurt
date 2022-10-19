#ifndef __XAVIER_COLLAB_GMIG_CORE_MANTEL_HPP__
#define __XAVIER_COLLAB_GMIG_CORE_MANTEL_HPP__

// STL
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
// xavier
#include "typedefs.hpp"

namespace xavier { namespace gmig { namespace core_mantle {
    
template<typename Scalar_, typename Index_=int>
struct Vertex {
    typedef Scalar_                         scalar_t;
    typedef Index_                          index_t;
    typedef nvis::fixed_vector<Scalar_, 3>  pos_t;
    typedef nvis::bounding_box<pos_t>       bbox_t;
    
    Vertex() : layer_id(-1), type(0) {}
    
    // helper functions
    scalar_t longitude() const { return position[1]; }
    scalar_t latitude() const { return position[0]; }
    scalar_t height() const { return position[2]; }
    index_t  sector_index() const { return long_sector*10+lat_sector; }
    bool valid() const { return fabs(type) == 1; }
    
    pos_t   position;
    index_t layer_id;
    index_t type;
    index_t lat_sector;
    index_t long_sector;
};

template<typename Scalar_, typename Index_ = int>
struct VertexDistance {
    typedef Scalar_                    scalar_t;
    typedef Index_                     index_t;
    typedef Vertex<scalar_t, index_t>  vertex_t;
    
    scalar_t operator()(const vertex_t& v1, const vertex_t& v2) const {
        return nvis::norm(v1.position - v2.position);
    }
};

template<typename Scalar_>
std::ostream& operator<<(std::ostream& os, const Vertex<Scalar_>& v) {
    os << '[' 
        << v.position    
        << ", layer="   << v.layer_id 
        << ", lat="     << v.latitude() 
        << ", long="    << v.longitude()
        << ", latsec="  << v.lat_sector 
        << ", longsec=" << v.long_sector
        << ", type="    << v.type 
        << ']';
    return os;
}

template<typename Scalar_, typename Index_ = int>
struct SectionID {
    typedef Index_                     index_t;
    typedef Scalar_                    scalar_t;
    typedef SectionID<Scalar_, Index_> self_type;
    typedef Vertex<Scalar_, Index_>    vertex_t;
    
    scalar_t longitude;
    index_t  lat_sector;
    index_t  long_sector;
    int      type;
    
    SectionID() : type(0) {} // invalid type to indicate non-initialized state
    
    SectionID(scalar_t l, index_t lats, index_t longs, int type_=1) 
        : longitude(l), lat_sector(lats), long_sector(longs), type(type_) {}
    
    SectionID(const vertex_t& v)
        : longitude(v.longitude()), lat_sector(v.lat_sector),
        long_sector(v.long_sector), type(v.type) {}
};

template<typename Scalar_, typename Index_ = int>  
inline bool operator<(const SectionID<Scalar_,Index_>& op1,
                      const SectionID<Scalar_,Index_>& op2) {
    if (op1.longitude < op2.longitude) return true;
    else if (op1.longitude > op2.longitude) return false;
    if (op1.long_sector < op2.long_sector) return true;
    else if (op1.long_sector > op2.long_sector) return false;
    if (op1.lat_sector < op2.lat_sector) return true;
    else if (op1.lat_sector > op2.lat_sector) return false;
    return op1.type < op2.type;
}

template<typename Scalar_, typename Index_ = int>  
inline bool operator==(const SectionID<Scalar_,Index_>& op1,
                       const SectionID<Scalar_,Index_>& op2)  {
    return op1.longitude == op2.longitude &&
           op1.long_sector == op2.long_sector &&
           op1.lat_sector == op2.lat_sector &&
           op1.type == op2.type;             
}

template<typename Scalar_, typename Index_ = int>
std::ostream& operator<<(std::ostream& os, 
                         const SectionID<Scalar_,Index_>& id) {
    os << "[long=" 
        << id.longitude << ", sec=[" 
        << id.lat_sector << ','
        << id.long_sector << "], type="
        << id.type
        << ']';
    return os;
}

template<typename Scalar_, typename Index_ = int>
struct LayerID : public SectionID<Scalar_, Index_> {
    typedef Index_                     index_t;
    typedef Scalar_                    scalar_t;
    typedef SectionID<Scalar_, Index_> base_type;
    typedef Vertex<Scalar_, Index_>    vertex_t;
    
    LayerID(const vertex_t& v)
        : base_type(v), layer_id(v.layer_id) {}
    
    int layer_id;
};

// indices of vertices contained on a given meridian
template<typename Scalar_, typename Index_ = int>
struct Section : public std::vector<Index_> {
    typedef Index_                     index_t;
    typedef std::vector<index_t>       base_t;
    typedef LayerID<Scalar_, Index_>   layer_id_t;
    
    Section() : base_t() {}
    Section(index_t id) : base_t(1, id) {}
    Section(size_t n, index_t id) : base_t(n, id) {}
    
    std::vector<layer_id_t> layers;
};


template<typename Scalar_, typename Index_ = int>  
inline bool operator<(const LayerID<Scalar_,Index_>& op1,
                      const LayerID<Scalar_,Index_>& op2) {
    typedef LayerID<Scalar_,Index_>        layer_id_t;
    typedef typename layer_id_t::base_type base_type;
    
    const base_type& sid1 = static_cast<const base_type&>(op1);
    const base_type& sid2 = static_cast<const base_type&>(op2);
    
    if (sid1<sid2) return true;
    else if (sid2<sid1) return false;
    else return op1.layer_id < op2.layer_id;
}

template<typename Scalar_, typename Index_ = int>  
inline bool operator==(const LayerID<Scalar_,Index_>& op1,
                      const LayerID<Scalar_,Index_>& op2)  {
    typedef LayerID<Scalar_,Index_>        layer_id_t;
    typedef typename layer_id_t::base_type base_type;
    
    // const base_type& sid1 = static_cast<base_type>(op1);
    // const base_type& sid2 = static_cast<base_type>(op2);
    
    return op1.base_type::operator==(op2) && op1.layer_id == op2.layer_id;             
}

template<typename Scalar_, typename Index_ = int>
std::ostream& operator<<(std::ostream& os, 
                         const LayerID<Scalar_,Index_>& id) {

    os << "[long=" 
        << id.longitude << ", sec=[" 
        << id.lat_sector << ','
        << id.long_sector << "], type="
        << id.type << ", layer_id="
        << id.layer_id 
        << ']';
    return os;
}

template<typename Scalar_, typename Index_ = int>
struct SectorID {
    typedef Index_                   index_t;
    typedef Scalar_                  scalar_t;
    typedef Vertex<Scalar_, Index_>  vertex_t;
    
    SectorID(const vertex_t& v) 
        : latitude_index(v.lat_sector), longitude_index(v.long_sector) {}
    
    index_t latitude_index;
    index_t longitude_index;
};

template<typename Scalar_, typename Index_ = int>
inline bool operator<(const SectorID<Scalar_, Index_>& op1,
                      const SectorID<Scalar_, Index_>& op2) {
    if (op1.longitude_index < op2.longitude_index) return true;
    else if (op1.longitude_index > op2.longitude_index) return false;
    else return op1.latitude_index < op2.latitude_index;
}

template<typename Scalar_, typename Index_ = int>
inline bool operator==(const SectorID<Scalar_, Index_>& op1,
                       const SectorID<Scalar_, Index_>& op2) {
    return op1.longitude_index == op2.longitude_index &&
           op1.latitude_index == op2.latitude_index;
}

template<typename Scalar_, typename Index_ = int>
std::ostream& operator<<(std::ostream& os, 
                         const SectorID<Scalar_, Index_>& s) {
    os << '[' << s.latitude_section << ", " << s.longitude_section << ']';          return os;
}

template<typename Scalar_, typename Index_ = int>
bool sanity_check(const std::vector<Vertex<Scalar_, Index_> >& 
                  vertices,
                  const Section<Scalar_, Index_>& sec, 
                  bool verbose=false) {
    typedef Scalar_ scalar_t;
    typedef Index_ index_t;
    typedef Section<Scalar_, Index_> section_t;
    typedef SectionID<Scalar_, Index_> section_id_t;
    typedef typename section_id_t::vertex_t vertex_t;
    
    std::set<scalar_t> longitudes;
    std::set<index_t> latitude_sec, longitude_sec;
    
    for (auto it=sec.begin(); it!= sec.end() ; ++it) {
        longitudes.insert(vertices[*it].longitude());
        latitude_sec.insert(vertices[*it].lat_sector);
        longitude_sec.insert(vertices[*it].long_sector);
    }
    
    // if (verbose) {
    //     std::cout << "there are " << longitudes.size()
    //         << " different longitudes in this section\n";
    //     std::cout << "there are " << latitude_sec.size()
    //         << " different latitude sections in this section\n";
    //     std::cout << "there are " << longitude_sec.size()
    //         << " different longitude sections in this section\n";
    // }
    
    return 
        longitudes.size() == 1 && 
        latitude_sec.size() == 1 && 
        longitude_sec.size() == 1;
}

// Mapping from Index_ to Set_<Value_>
template<typename Index_, typename Set_, 
         typename Value_ = typename Set_::value_type>
class SetContainer : public std::map<Index_, Set_> {
public:
    typedef Index_ index_type;
    typedef Value_ value_type;
    typedef Set_   set_type;
    typedef typename std::map<Index_, Set_> base_type;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    
public:
    SetContainer() : base_type() {}
    
    value_type& back(const index_type& id) {
        return base_type::operator[](id).back();
    }
    
    value_type& front(const index_type& id) {
        return base_type::operator[](id).front();
    }
    
    const value_type& back(const index_type& id) const {
        return base_type::operator[](id).back();
    }
    
    const value_type& front(const index_type& id) const {
        return base_type::operator[](id).front();
    }
    
    // overload map::insert function with different signature
    using base_type::insert; // prevent name hiding
    bool insert(const index_type& id, const value_type& v) {
        // we attempt to create a new set with index id
        // that we initialize with a single element v
        std::pair<iterator, bool> where = 
            base_type::insert(std::make_pair(id, set_type(1, v)));
        if (!where.second) {
            // insertion failed because a set with index id already exists
            where.first->second.push_back(v);
            return false;
        }
        return true;
    }
    
    bool exists(const index_type& id) const {
        return base_type::find(id) != base_type::end();
    }
};



} // namespace core_mantle
} // namespace gmig
} // namespace xavier

#endif
