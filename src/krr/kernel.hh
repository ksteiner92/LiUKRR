/*
 * kernel.hh
 *
 *  Created on: 30.10.2016
 *      Author: Klaus Steiner
 */

#ifndef KERNEL_HH_
#define KERNEL_HH_

#include "types.h"

template<class T> struct KernelType {
    typedef std::function<
            double(const T* desc1, const T* desc2, const size_t n,
                    const double sigma)> type;
};

template<class T> struct Kernel {

    /**
     * Laplace kernel with L1 norm
     */
    static const constexpr auto LAPLACE =
            STATIC_LAMBDA(const T* desc1, const T* desc2, const size_t n, const double sigma) {
                T* diff = new T[n * n];
                MKLvSub<T>::call(n * n, desc1, desc2, diff);
                double norm1 = BlasAsum<T>::call(n * n, diff, 1);
                delete[] diff;
                return exp(- norm1 / sigma);
            };
};

#endif /* KERNEL_HH_ */
