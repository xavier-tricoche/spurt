
#ifndef INCLUDE_KDTREE_NODE2_HPP
#define INCLUDE_KDTREE_NODE2_HPP

#ifdef KDTREE_DEFINE_OSTREAM_OPERATORS
#include <ostream>
#endif

#include <cstddef>
#include <cmath>
#include <vector>
#include <data/heap.hpp>

#include "node.hpp"

namespace KDTree
{

    /*! Find the N nearest nodes to __val from __node

      If many nodes are equidistant to __val, the node with the lowest memory
      address is returned.

      \return the nearest node of __end node if no nearest node was found for the
      given arguments.
     */
    template <class SearchVal,
              typename _Val, typename _Cmp,
              typename _Acc, typename _Dist,
              typename _Predicate>
    inline std::vector< std::pair< const _Node<_Val> *,
                                   std::pair< size_t, 
                                              typename _Dist::distance_type
                                            >
                                 >
                      >
    _S_node_n_nearest(const size_t __k, const size_t __n, 
                      size_t __dim, SearchVal const &__val,
                      bool node_is_candidate,
                      const _Node<_Val> *__node, const _Node_base *__end,
                      const _Node<_Val> *__best, 
                      typename _Dist::distance_type __max,
                      const _Cmp &__cmp, const _Acc &__acc, const _Dist &__dist,
                      _Predicate __p)
    {
        const _Node_base *pcur = __node;
        const _Node_base *cur = _S_node_descend(__dim % __k, __cmp, __acc, 
                                                __val, __node);
        size_t cur_dim = __dim + 1;

        typedef _Node<_Val> node_type;
        typedef const node_type* const_node_ptr; 
        typedef std::pair< size_t, typename _Dist::distance_type > dim_dist;
        typedef std::pair< const_node_ptr, dim_dist > nn_type;

        struct Closer {
            bool operator()(const nn_type& n0, const nn_type& n1) const {
                return n0.second.second < n1.second.second;
            }
        };
        Closer less_;

        spurt::bounded_heap<nn_type, Closer> n_nn(__n);

        if (node_is_candidate) 
            n_nn.push(nn_type(__best, dim_dist(__dim, __max)));

        // find the n smallest distances <= __max that satisfy __p 
        // in direct descent
        while (cur)
        {
            if (__p(static_cast<const_node_ptr>(cur)->_M_value))
            {
                typename _Dist::distance_type d = 0;
                for (size_t i = 0; i != __k; ++i)
                    d += _S_node_distance(i, __dist, __acc, __val, static_cast<const_node_ptr>(cur)->_M_value);
                d = sqrt(d);

                nn_type candidate(static_cast<const_node_ptr>(cur), 
                                  dim_dist(cur_dim, d));
                if (d <= __max)
                {
                    // will trim top is capacity exceeded
                    n_nn.push(candidate, false); 
                    // if we already have n neighbors, largest neighbor
                    // distance is new max
                    if (n_nn.full()) __max = n_nn.top().second.second;
                }
            }
            pcur = cur;
            cur = _S_node_descend(cur_dim % __k, __cmp, __acc, __val, cur);
            ++cur_dim;
        }
        // Swap cur to prev, only prev is a valid node.
        cur = pcur;
        --cur_dim;
        pcur = NULL;
        // Probe all node's children not visited yet (siblings of the visited nodes).
        const _Node_base *probe = cur;
        const _Node_base *pprobe = probe;
        const _Node_base *near_node;
        const _Node_base *far_node;
        size_t probe_dim = cur_dim;
        if (_S_node_compare(probe_dim % __k, __cmp, __acc, __val,       
                            static_cast<const_node_ptr>(probe)->_M_value))
            near_node = probe->_M_right;
        else
            near_node = probe->_M_left;
        if (near_node
            // only visit node's children if node's plane intersect hypersphere
            && (sqrt(_S_node_distance(probe_dim % __k, __dist, __acc, __val, 
                                      static_cast<const_node_ptr>(probe)
                                      ->_M_value)) <= __max))
        {
            probe = near_node;
            ++probe_dim;
        }
        while (cur != __end)
        {
            while (probe != cur)
            {
                if (_S_node_compare(probe_dim % __k, __cmp, __acc, __val, 
                                    static_cast<const_node_ptr>(probe)
                                    ->_M_value))
                {
                    near_node = probe->_M_left;
                    far_node = probe->_M_right;
                }
                else
                {
                    near_node = probe->_M_right;
                    far_node = probe->_M_left;
                }
                if (pprobe == probe->_M_parent) // going downward ...
                {
                    if (__p(static_cast<const_node_ptr>(probe)->_M_value))
                    {
                        typename _Dist::distance_type d = 0;
                        for (size_t i = 0; i < __k; ++i)
                            d += _S_node_distance(i, __dist, __acc, __val, 
                                                  static_cast<const_node_ptr>
                                                  (probe)->_M_value);
                        d = sqrt(d);
                        nn_type new_candidate(
                            static_cast<const_node_ptr>(probe), 
                            dim_dist(probe_dim, d));
                        if (d < __max)
                        {
                            // will trim top is capacity exceeded
                            n_nn.push(new_candidate, false);
                            // if we already have n neighbors, largest neighbor
                            // distance is new max
                            if (n_nn.full()) __max = n_nn.top().second.second;
                        }
                    }
                    pprobe = probe;
                    if (near_node)
                    {
                        probe = near_node;
                        ++probe_dim;
                    }
                    else if (far_node &&
                             // only visit node's children if node's plane intersect hypersphere
                             sqrt(_S_node_distance(probe_dim % __k, __dist, __acc, __val, static_cast<const _Node<_Val> *>(probe)->_M_value)) <= __max)
                    {
                        probe = far_node;
                        ++probe_dim;
                    }
                    else
                    {
                        probe = probe->_M_parent;
                        --probe_dim;
                    }
                }
                else // ... and going upward.
                {
                    if (pprobe == near_node && far_node
                        // only visit node's children if node's plane intersect hypersphere
                        && sqrt(_S_node_distance(probe_dim % __k, __dist, 
                                                 __acc, __val, 
                                                 static_cast<const_node_ptr>(probe)->_M_value)) <= __max)
                    {
                        pprobe = probe;
                        probe = far_node;
                        ++probe_dim;
                    }
                    else
                    {
                        pprobe = probe;
                        probe = probe->_M_parent;
                        --probe_dim;
                    }
                }
            }
            pcur = cur;
            cur = cur->_M_parent;
            --cur_dim;
            pprobe = cur;
            probe = cur;
            probe_dim = cur_dim;
            if (cur != __end)
            {
                if (pcur == cur->_M_left)
                    near_node = cur->_M_right;
                else
                    near_node = cur->_M_left;
                if (near_node
                    // only visit node's children if node's plane intersect hypersphere
                    && (sqrt(_S_node_distance(cur_dim % __k, __dist, __acc, __val, static_cast<const _Node<_Val> *>(cur)->_M_value)) <= __max))
                {
                    probe = near_node;
                    ++probe_dim;
                }
            }
        }
        return n_nn;
    }

} // namespace KDTree

#endif