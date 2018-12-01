/*
 * machine.cc
 *
 *  Created on: 18.09.2016
 *      Author: Klaus Steiner
 */

#include "machine.h"

template<class T> Machine<T>::Machine(Descriptor<T> *descriptor,
        const typename KernelType<T>::type &kernel) :
        descriptor(descriptor), alphas(NULL), sigma(4.0e4), lambda(1.0e-4), includeElements(
        NULL), numAtoms(0), kernelFunc(kernel) {
}

template<class T> Machine<T>::~Machine() {
    if (alphas != NULL)
        delete[] alphas;
    if (includeElements != NULL)
        delete[] includeElements;
}

template<class T> void Machine<T>::setSigma(const double sigma) {
    this->sigma = sigma;
}

template<class T> const double Machine<T>::getSigma() const {
    return sigma;
}

template<class T> void Machine<T>::setLambda(const double lambda) {
    this->lambda = lambda;
}

template<class T> const double Machine<T>::getLambda() const {
    return lambda;
}

template<class T> void Machine<T>::train(DataSet<T> *dataSet) {
    matrices.clear();
    descriptor->calculate(dataSet, matrices);
    /* calculate the kernel */
    const auto begin = std::chrono::steady_clock::now();
    std::cout << "[Info] Calculating kernel... ";
    std::cout.flush();
    lapack_int n = matrices.size();
    numAtoms = dataSet->getNumAtomsWithRestriction();
    double* kernel = new double[n * n];
    const size_t kernelDims[2] = { matrices.size(), matrices.size() };
    size_t kernelPos[2];
    double val;
//#pragma omp parallel for shared(val)
    for (size_t i = 0; i < n; i++) {
        kernelPos[0] = i;
        kernelPos[1] = kernelPos[0];
        kernel[getArrayPosRowMajor(kernelPos, kernelDims, 2)] = 1.0 + lambda;
        for (kernelPos[1] = kernelPos[0] + 1; kernelPos[1] < n;
                kernelPos[1]++) {
            val = kernelFunc(matrices[kernelPos[0]].m.get(),
                    matrices[kernelPos[1]].m.get(), numAtoms, sigma);
//#pragma omp atomic write
            kernel[kernelPos[0] * n + kernelPos[1]] = val;
        }
    }
    const auto end = std::chrono::steady_clock::now();
    std::cout << "[OK] Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                    end - begin).count() << "ms" << std::endl;
    /* Solve the equation system */
    std::cout << "[Info] Calculating alphas... ";
    std::cout.flush();
    lapack_int ipiv[n];
    lapack_int info = 0;
    if (alphas != NULL)
        delete[] alphas;
    alphas = new double[n];
    std::fill_n(alphas, n, 0.0);
    T *tmp = new T[n];
    dataSet->getEnergies(tmp);
    double *energies = ConvertToDoubleArray<T>::convert(tmp, n);
    info = LapackSsytrf<double>::call(LAPACK_ROW_MAJOR, 'U', n, kernel, n,
            ipiv);
    if (info < 0) {
        std::stringstream ss;
        ss << "Matrix factorization: " << "the " << -info
                << "-th parameter had an illegal value.";
        throw ss.str();
    } else if (info > 0) {
        throw "Matrix factorization: Matrix singular.";
    }
    info = LapackSytri<double>::call(LAPACK_ROW_MAJOR, 'U', n, kernel, n, ipiv);
    if (info < 0) {
        std::stringstream ss;
        ss << "the " << -info << "-th parameter had an illegal value.";
        throw ss.str();
    } else if (info > 0) {
        throw "Matrix singular.";
    }
    BlasSsymv<double>::call(CblasRowMajor, CblasUpper, n, 1.0, kernel, n,
            energies, 1, 0.0, alphas, 1);
    const auto end2 = std::chrono::steady_clock::now();
    std::cout << "[OK] Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count()
            << "ms" << std::endl;
    std::cout.flush();
    delete[] energies;
    delete[] kernel;
}

template<class T> void Machine<T>::predict(DataSet<T> *dataSet,
        PredictionResults &results) const {
    std::cout << "[Info] Predicting energies... ";
    std::cout.flush();
    const auto begin = std::chrono::steady_clock::now();
    std::vector<DescriptorMatrix<T>> predDescs;
    descriptor->calculate(dataSet, predDescs);
    T *energies = new T[predDescs.size()];
    const size_t numAtoms = dataSet->getNumAtomsWithRestriction();
    dataSet->getEnergies(energies);
    results.mae = 0.0;
    results.maePerAtom = 0.0;
    for (size_t i = 0; i < predDescs.size(); i++) {
        Prediction pred;
        pred.name = predDescs[i].name;
        pred.numAtoms = predDescs[i].numAtoms;
        pred.value = 0.0;
        for (size_t j = 0; j < matrices.size(); j++) {
            pred.value += alphas[j]
                    * kernelFunc(predDescs[i].m.get(), matrices[j].m.get(),
                            numAtoms, sigma);
        }
        pred.valuePerAtom = pred.value / predDescs[i].numAtoms;
        pred.error = energies[i] - pred.value;
        pred.errorPerAtom = pred.error / predDescs[i].numAtoms;
        results.mae += fabs(pred.error);
        results.maePerAtom += fabs(pred.errorPerAtom);
        results.results.push_back(pred);
    }
    delete[] energies;
    results.mae /= predDescs.size();
    results.maePerAtom /= predDescs.size();
    const auto end = std::chrono::steady_clock::now();
    std::cout << "[OK] Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                    end - begin).count() << "ms" << std::endl;
    std::cout.flush();
}

/**
 * TODO save the state of the machine
 * save file must contain:
 * - alphas
 * - sigma
 * - lambda
 * - descriptors on which was trained
 * - descriptor type
 * - kernel type / kernel function
 * save file could contain?
 * - sub set informations
 */
template<class T> void Machine<T>::save(const std::string &fileName) const {
    if (alphas == NULL)
        return;
    /* Write out the solution in a new mat file */
    std::cout << "[OK]" << std::endl << "[Info] Writing solution in \""
            << fileName << "\"... ";
    MatFileHandle* matfp = Mat_CreateVer(fileName.c_str(), NULL, MAT_FT_MAT5);
    if (matfp == NULL) {
        throw "Could not create output file!";
    }
    const size_t numCompounds = matrices.size();
    size_t dimsAlphas[2] = { 1, numCompounds };
    MatVariable *matvar = Mat_VarCreate("A", MAT_C_SINGLE, MAT_T_SINGLE, 2,
            dimsAlphas, alphas, 0);
    if (matvar == NULL) {
        throw "Could not create 'A' variable in output file!";
    }
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);

//MatVariable *matlabel = Mat_VarCreate("desc", MAT_C_CHAR, MAT_T_STRING, )
    size_t dimC[3] = { numCompounds, numAtoms, numAtoms };
    T* x = new T[numCompounds * numAtoms * numAtoms];
    matvar = Mat_VarCreate("C", MAT_C_SINGLE, MAT_T_SINGLE, 3, dimC, x, 0);
    if (matvar == NULL) {
        delete[] x;
        throw "Could not create 'A' variable in output file!";
    }
    size_t pos[3];
    size_t pos2[2];
    size_t dimC2[2] = { numAtoms, numAtoms };
    pos[0] = 0;
    for (size_t i = 0; i < numCompounds; i++) {
        for (pos[1] = 0; pos[1] < numAtoms; pos[1]++) {
            for (pos[2] = 0; pos[2] < numAtoms; pos[2]++) {
                pos2[0] = pos[1];
                pos2[1] = pos[2];
                /*x[getArrayPosRowMajor(pos, dimC, 3)] =
                 (pos[1] == pos[2] ?
                 desc.diag[pos[1]] :
                 desc.offDiag[getArrayPosRowMajor(pos2, dimC2, 2)]);*/
            }
        }
        pos[0]++;
    }
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    std::cout << "[OK]" << std::endl;
}

template<class T> void Machine<T>::load(const std::string &fileName) {

}

template class Machine<float> ;
template class Machine<double> ;
