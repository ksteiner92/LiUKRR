/*
 * descriptor.h
 *
 *  Created on: 19.09.2016
 *      Author: Klaus Steiner
 */

#ifndef DESCRIPTOR_H_
#define DESCRIPTOR_H_

#include <iostream>
#include <vector>
#include <typeinfo>
#include <memory>
#include <functional>

#include "dataset.h"
#include "descriptor.h"
#include "types.h"
#include "sorting.hh"

template<class T> class DataSet;

template<class T> class Descriptor {
public:
    Descriptor();

    virtual ~Descriptor() {
    }

    /**
     * Sets the sorting of the resulting descriptor matrix. The algorithm itself
     * is provided as a lambda function.
     *
     * @param sorting The sorting algorithm as lambda function.
     */
    void setSorting(const typename SortingAlgorithmType<T>::type &sorting);

    /**
     * Calculates a descriptor matrices for a given data set and stores the results
     * in a vector. For every matrix the sorting algorithm is executed. Default is
     * no sorting.
     *
     * @param dataSet The data (sub) set for which the descriptor matrices should be
     * calculated.
     * @param descriptors A vector in which the results are stored
     */
    virtual void calculate(DataSet<T> *dataSet,
            std::vector<DescriptorMatrix<T>> &descriptors) = 0;

protected:
    typename SortingAlgorithmType<T>::type sorting;
};

template<class T> class CoulombDescriptor: public Descriptor<T> {
public:
    void calculate(DataSet<T> *dataSet,
            std::vector<DescriptorMatrix<T>> &descriptors);
};

template<class T> class SineDescriptor: public Descriptor<T> {
public:
    void calculate(DataSet<T> *dataSet,
            std::vector<DescriptorMatrix<T>> &descriptors);
};

#endif /* DESCRIPTOR_H_ */
