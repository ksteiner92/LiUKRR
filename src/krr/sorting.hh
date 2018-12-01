/*
 * sorting.hh
 *
 *  Created on: 30.10.2016
 *      Author: Klaus Steiner
 */

#ifndef SORTING_HH_
#define SORTING_HH_

#include "types.h"

template<class T> struct SortingAlgorithmType {
    typedef std::function<void(T* &descriptor, const size_t n)> type;
};

template<class T> struct SortingAlgorithm {
    static const constexpr auto ROW_L1_DESC =
            STATIC_LAMBDA(T* &descriptor, const size_t n) {
                std::vector<size_t> row_order;
                for(size_t i = 0; i < n; i++)
                    row_order.push_back(i);
                T asum[n];
                std::fill_n(asum, n, -1.0);
                std::sort(row_order.begin(), row_order.end(), [n, descriptor, &asum](size_t a, size_t b) {
                            return (asum[a] >= 0 ? asum[a] : (asum[a] = BlasAsum<T>::call(n, descriptor + a * n, 1))) >
                            (asum[b] >= 0 ? asum[b] : (asum[b] = BlasAsum<T>::call(n, descriptor + b * n, 1)));});
                T* tmp = new T[n * n];
                for(size_t i = 0; i < n; i++) {
                    size_t idx = row_order[i];
                    for (size_t j = i; j < n; j++) {
                        tmp[i * n + j] = descriptor[row_order[j] * n + idx];
                        tmp[j * n + i] = tmp[i * n + j];
                    }
                }
                std::copy(tmp, tmp + n * n, descriptor);
                delete[] tmp;
            };

    static const constexpr auto ROW_L2_DESC =
            STATIC_LAMBDA(T* &descriptor, const size_t n) {
                std::vector<size_t> row_order;
                for(size_t i = 0; i < n; i++)
                    row_order.push_back(i);
                T norm_2[n];
                std::fill_n(norm_2, n, -1.0);
                std::sort(row_order.begin(), row_order.end(), [n, descriptor, &norm_2](size_t a, size_t b) {
                            return (norm_2[a] >= 0 ? norm_2[a] : (norm_2[a] = BlasNrm2<T>::call(n, descriptor + a * n, 1) )) >
                            (norm_2[b] >= 0 ? norm_2[b] : (norm_2[b] = BlasNrm2<T>::call(n, descriptor + b * n, 1)));});
                T* tmp = new T[n * n];
                for(size_t i = 0; i < n; i++) {
                    size_t idx = row_order[i];
                    for (size_t j = i; j < n; j++) {
                        tmp[i * n + j] = descriptor[row_order[j] * n + idx];
                        tmp[j * n + i] = tmp[i * n + j];
                    }
                }
                std::copy(tmp, tmp + n * n, descriptor);
                delete[] tmp;
            };
};

#endif /* SORTING_HH_ */
