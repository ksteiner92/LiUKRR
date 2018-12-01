/*
 * main.cc
 *
 *  Created on: 06.09.2016
 *      Author: Klaus Steiner
 */

#include <iostream>
#include <algorithm>
#include <limits.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <map>

#include "progargparser.h"
#include "dataset.h"
#include "descriptor.h"
#include "machine.h"
#include "types.h"
#include "quantity.h"

void printHelp(const std::string &programName) {
    std::cout << "Usage: " << programName
            << " <Program Mode> [Options] <Input> <Output>" << std::endl;
    std::cout << "Program Modes:" << std::endl;
    std::cout
            << "  --train:\tTrains the program for a given data set in <Input> and writes the trained machine in <Output>"
            << std::endl;
    std::cout
            << "  --compute:\tComputes The formation energy for a given data set in <Input>"
            << std::endl;
}

int main(int argc, char **argv) {
    /*QMDataSet_7 *data = new QMDataSet_7(
            "/home/klaus/dev/LiU/data/qm7.mat");
    data->setSubSet(2000,  2001);
    Descriptor<float> *descriptor = new CoulombDescriptor<float>;

    descriptor->setSorting(SortingAlgorithm<float>::ROW_L1_DESC);
    Machine<double> machine(descriptor, Kernel<float>::LAPLACE);
    machine.setSigma(4.0e4);
    machine.setLambda(1.0e-4);
    machine.train(data);
    data->setSubSetToAllExcludedElements();*/

    //VASPDataSet *data = new VASPDataSet("/home/klaus/dev/LiU/data/mp-basicStructures_Symmetrized", true);
    //DataSet<float> *data = new QMDataSet_7(inputFile);
    //Descriptor<double> *descriptor = new SineDescriptor<double>;
    //Descriptor<float> *descriptor = new CoulombDescriptor<float>;
    //descriptor->setSorting(SortingAlgorithm<double>::ROW_L1_DESC);
    //data->setSubSet(1000,  1001);
    //data->setSubSet(0,  3000);
    /*Structures<float> structures;
     data->parseStructures(structures);
     std::vector<DescriptorMatrix<float>> descriptors;
     descriptor->calculate(data, descriptors);
     for (int i = 0; i < 2; i++) {
     std::cout << i << ": num atoms struct: " << structures.structures[i].atoms.size() << std::endl;
     std::cout << i << ": num atoms desc:   " << descriptors[i].numAtoms << std::endl;
     }
     VASPDataSet *data2 = new VASPDataSet("/home/klaus/dev/LiU/data/mp-basicStructures_Symmetrized");
     data2->setSubSet(0, 10);
     Structures<double> structures2;
     data2->parseStructures(structures2);
     for (int i = 0; i < 10; i++)
     std::cout << i << ": num atoms struct: " << structures2.structures[i].atoms.size() << std::endl;*/

    /*Machine<double> machine(descriptor, Kernel<double>::LAPLACE);
     machine.setSigma(4.0e4);
     machine.setLambda(1.0e-4);
     machine.train(data);
     data->setSubSetToAllExcludedElements();
     PredictionResults results;
     machine.predict(data, results);
     std::cout << "\tMAE: " << results.mae << " kcal/mol" << std::endl;
     std::cout << "\tMAE: " << results.mae * 0.043e3 << " meV"
     << std::endl;
     std::cout << "\tMAE per atom: " << results.maePerAtom << " kcal/mol"
     << std::endl;
     std::cout << "\tMAE per atom: " << results.maePerAtom * 0.043e3
     << " meV" << std::endl;
     delete data;
     delete descriptor;*/

    //error(inputFile);
    /*const auto begin = std::chrono::steady_clock::now();
     QMDataSet_7 *data = new QMDataSet_7(inputFile);
     //data->setSubSet(0, 7000);
     CoulombDescriptor<float> *descriptor = new CoulombDescriptor<float>;
     QuantumMachine<float> qm(descriptor);
     //data->setSubSet(0, 7000);
     data->setRandomSubSet(5000);
     qm.train(data);*/
    /*const QuantumMachine<float>::KernelFunctionType gaussKernel =
     [](const float* desc1, const float* desc2, const size_t n, const double sigma) {
     float* diff = new float[n * n];
     MKLvSub<float>::call(n * n, desc1, desc2, diff);
     float* work = new float[n];
     MKL_INT mkl_n = n;
     double norm2 = LapackLansy<float>::call(LAPACK_ROW_MAJOR, 'M', 'U', mkl_n, diff, mkl_n);
     delete[] work;
     delete[] diff;
     return exp(- norm2 * norm2 / (2 * sigma * sigma));
     };*/
    //qm.setKernel(gaussKernel);
    //data->setSubSet(7000, 7165);
    /*data->setSubSetToAllExcludedElements();
     std::vector<float> results;
     qm.predict(data, results);
     //qm.save(outputFile);
     delete descriptor;
     delete data;
     const auto end = std::chrono::steady_clock::now();
     std::cout << "Time: " << std::chrono::duration_cast
     < std::chrono::seconds
     > (end - begin).count() << "s" << std::endl;*/
    /*char result[ PATH_MAX];
     ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
     VASPDataSet *vaspData = new VASPDataSet(
     std::string(result, (count > 0) ? count : 0));
     vaspData->preProcess(inputFile);*/
    return 0;
}
