/*
 * machine.h
 *
 *  Created on: 18.09.2016
 *      Author: Klaus Steiner
 */

#ifndef MACHINE_H_
#define MACHINE_H_

#include <algorithm>
#include <iostream>
#include <chrono>
#include <math.h>
#include <memory>
#include <sstream>
#include <functional>
#include <iomanip>
#include <fstream>

#include "descriptor.h"
#include "dataset.h"
#include "types.h"
#include "kernel.hh"

template<class T> class DataSet;
template<class T> class Descriptor;

struct Prediction {
    double value;
    double valuePerAtom;
    double error;
    double errorPerAtom;
    size_t numAtoms;
    std::string name;
};

struct PredictionResults {
    double mae;
    double maePerAtom;
    std::vector<Prediction> results;
};

template<class T> class Machine {
public:
    /**
     *
     */
    Machine(Descriptor<T> *descriptor,
            const typename KernelType<T>::type &kernel);

    /**
     *
     */
    ~Machine();

    void setSigma(const double sigma);

    const double getSigma() const;

    void setLambda(const double lambda);

    const double getLambda() const;

    /**
     *
     */
    void train(DataSet<T> *dataSet);

    /**
     *
     */
    void predict(DataSet<T> *dataSet, PredictionResults &results) const;

    /**
     *
     */
    void save(const std::string &fileName) const;

    /**
     *
     */
    void load(const std::string &fileName);

private:
    Descriptor<T> *descriptor;
    std::vector<DescriptorMatrix<T>> matrices;
    size_t numAtoms;
    typename KernelType<T>::type kernelFunc;
    double *alphas;
    bool* includeElements;
    double sigma;
    double lambda;

};

#endif /* MACHINE_H_ */
