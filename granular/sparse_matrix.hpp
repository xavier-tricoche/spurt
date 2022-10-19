#include <set>
#include <list>
#include <vector>

// #define DEBUG

template<typename T>
struct sparse_matrix {

    sparse_matrix() {}
    
    void initialize(const std::vector<std::list<std::pair<unsigned int, T> > >& rows) {
        typedef std::pair<unsigned int, T>          pair_type;
        typedef std::list<pair_type>                list_type;
        typedef typename list_type::const_iterator  iterator_type;
        
        __rows.resize(rows.size() + 1);
        
        int n = 0;
        for (int i = 0 ; i < rows.size() ; ++i) {
            n += rows[i].size() + 1;
        }
        __cols.resize(n);
        
        int offset = 0;
        for (int i = 0 ; i < rows.size() ; ++i) {
            __rows[i] = offset;
            const list_type& column = rows[i];
            for (iterator_type it = column.begin() ; it != column.end() ; ++it) {
                __cols[offset++] = *it;
            }
        }
        __rows.back() = __cols.size();
    }
    
    T operator()(int i, int j) const {
    
#ifdef DEBUG
        std::cerr << "calling operator()(" << i << ", " << j << ")\n";
#endif
        
        unsigned int row_starts = __rows[i];
        unsigned int row_ends   = __rows[i+1];
        // check if row is empty
        if (row_starts == row_ends) {
            return T(0);
        }
        
#ifdef DEBUG
        std::cerr << "row starts at " << row_starts << " and ends at " << row_ends << '\n';
#endif
        
        unsigned int first_col = __cols[row_starts].first;
        unsigned int last_col = __cols[row_ends-1].first;
        
#ifdef DEBUG
        std::cerr << "column ids range from " << first_col << " to " << last_col << '\n';
#endif
        
        // check if column is out of range
        if (j < first_col || j > last_col) {
            return T(0);
        } else if (j == first_col) {
            return __cols[row_starts].second;
        } else if (j == last_col) {
            return __cols[row_ends-1].second;
        }
        
        // start binary search
        unsigned int min_id = row_starts;
        unsigned int max_id = row_ends - 1;
        while (min_id + 1 < max_id) {
#ifdef DEBUG
            std::cerr << "min=" << min_id << ", max=" << max_id << '\n';
#endif
            unsigned int mid_id = (min_id + max_id) / 2;
            unsigned int col = __cols[mid_id].first;
            if (col == j) {
                return __cols[mid_id].second;
            } else if (col > j) {
                max_id = mid_id;
            } else {
                min_id = mid_id;
            }
        }
        
        return T(0);
    }
    
    void row_entries(std::vector<std::pair<unsigned int, T> >& row, unsigned int r) const {
        row.clear();
        for (unsigned int i = __rows[r] ; i != __rows[r+1] ; ++i) {
            row.push_back(__cols[i]);
        }
    }
    
    std::vector<std::pair<unsigned int, T> >    __cols;
    std::vector<unsigned int>                   __rows;
    
};

