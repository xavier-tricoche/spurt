#include <data/heap.hpp>
#include <iostream>


template<typename T, typename Comp> 
std::ostream& operator<<(std::ostream& os, const spurt::heap<T, Comp>& h) {
    std::copy(h.begin(), h.end(), std::ostream_iterator<int>(os, ", "));
    return os;
}


int main(int argc, char* argv[]) {
    int values[15] = { 0, 4, 2, 7, 38, 39, 21, 24, 56, 102, 324, 54, 23, 12, 1025 };

    typedef spurt::heap<int> max_heap_t;
    typedef spurt::bounded_heap<int> max_bheap_t;

    struct invless {
        bool operator()(int a, int b) {
            return b < a;
        }
    };

    invless _less;

    typedef spurt::heap<int, invless> min_heap_t;
    typedef spurt::bounded_heap<int, invless> min_bheap_t;

    max_heap_t max_heap(values, values+11);
    min_heap_t min_heap(values, values+11);
    max_bheap_t max_bheap(values, values+10, 10);
    min_bheap_t min_bheap(values, values+10, 10);

    std::cout << "max_heap is\n" << max_heap << '\n';
    std::cout << "is_heap: " << std::is_heap(max_heap.begin(), max_heap.end()) << '\n';
    std::cout << "min heap is\n"
              << min_heap << '\n';
    std::cout << "is_heap: " << std::is_heap(min_heap.begin(), min_heap.end(), _less) << '\n';
    std::cout << "bounded max heap is\n"
              << max_bheap << '\n';
    std::cout << "is_heap: " << std::is_heap(max_bheap.begin(), max_bheap.end()) << '\n';
    std::cout << "bounded min heap is\n"
              << min_bheap << '\n';
    std::cout << "is_heap: " << std::is_heap(min_bheap.begin(), min_bheap.end(), _less) << '\n';

    std::cout << "popping max heap returns " << max_heap.pop() << '\n';
    std::cout << "max heap is now\n" << max_heap << '\n';

    std::cout << "is_heap: " << std::is_heap(max_heap.begin(), max_heap.end()) << '\n';

    std::cout << "popping min heap returns " << min_heap.pop() << '\n';
    std::cout << "min heap is now\n"
              << min_heap << '\n';
    std::cout << "is_heap: " << std::is_heap(min_heap.begin(), min_heap.end(), _less) << '\n';

    std::cout << "pushing 253 on max heap gives\n";
    max_heap.push(253);
    std::cout << max_heap << '\n';
    std::cout << "is_heap: " << std::is_heap(max_heap.begin(), max_heap.end()) << '\n';

    std::cout << "pushing 253 on min heap gives\n";
    min_heap.push(253);
    std::cout << min_heap << '\n';
    std::cout << "is_heap: " << std::is_heap(min_heap.begin(), min_heap.end(), _less) << '\n';

    std::cout << "pushing 253 on max bounded heap and trimming bottom gives\n";
    max_bheap.push(253);
    std::cout << max_bheap << '\n';
    std::cout << "is_heap: " << std::is_heap(max_bheap.begin(), max_bheap.end()) << '\n';
    std::cout << "pushing 17 on max bounded heap and trimming top gives\n";
    max_bheap.push(17, false);
    std::cout << max_bheap << '\n';
    std::cout << "is_heap: " << std::is_heap(max_bheap.begin(), max_bheap.end()) << '\n';

    std::cout << "pushing 253 on min bounded heap and trimming bottom gives\n";
    min_bheap.push(253);
    std::cout << min_bheap << '\n';
    std::cout << "is_heap: " << std::is_heap(min_bheap.begin(), min_bheap.end(), _less) << '\n';
    std::cout << "pushing 17 on min bounded heap and trimming top gives\n";
    min_bheap.push(17, false);
    std::cout << min_bheap << '\n';
    std::cout << "is_heap: " << std::is_heap(min_bheap.begin(), min_bheap.end(), _less) << '\n';
}