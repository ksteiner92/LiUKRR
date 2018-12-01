/*
 * descriptor.cc
 *
 *  Created on: 19.09.2016
 *      Author: Klaus Steiner
 */

#include "descriptor.h"

template<class T> Descriptor<T>::Descriptor() :
        sorting([](T* &descriptor, const size_t n) {}) {
}

template<class T> void Descriptor<T>::setSorting(
        const typename SortingAlgorithmType<T>::type &sorting) {
    this->sorting = sorting;
}

template<class T> void CoulombDescriptor<T>::calculate(DataSet<T> *dataSet,
        std::vector<DescriptorMatrix<T>> &descriptors) {
    QMDataSet_7 *qm7DataSet = dynamic_cast<QMDataSet_7*>(dataSet);
    if (qm7DataSet != NULL) {
        qm7DataSet->parseCoulombMatrices(descriptors);
    } else {
        Structures<T> structures;
        dataSet->getStructures(structures);
        const size_t numAtoms = dataSet->getNumAtomsWithRestriction();
        for (const Structure<T> structure : structures.structures) {
            T* matrix = new T[numAtoms * numAtoms];
            std::fill_n(matrix, numAtoms * numAtoms, 0.0);
            DescriptorMatrix<T> descriptor;
            descriptor.numAtoms = structure.atoms.size();
            for (size_t i = 0; i < descriptor.numAtoms; i++) {
                if (structure.atoms[i].Z > 0) {
                    matrix[i * numAtoms + i] = 0.5
                            * pow(structure.atoms[i].Z, 2.4);
                    for (size_t j = i + 1; j < descriptor.numAtoms; j++) {
                        T phi = 0.0;
                        for (int k = 0; k < 3; k++) {
                            const T diff = structure.atoms[i].r[k]
                                    - structure.atoms[j].r[k];
                            phi += diff * diff;
                        }
                        phi = 1.0 / sqrt(phi);
                        matrix[i * numAtoms + j] = structure.atoms[i].Z
                                * structure.atoms[j].Z * phi;
                        matrix[j * numAtoms + i] = matrix[i * numAtoms + j];
                    }
                }
            }
            this->sorting(matrix, numAtoms);
            descriptor.name = structure.name;
            descriptor.m = std::shared_ptr<T>(matrix);
            descriptors.push_back(descriptor);
        }
    }
}

template<class T> void SineDescriptor<T>::calculate(DataSet<T> *dataSet,
        std::vector<DescriptorMatrix<T>> &descriptors) {
    Structures<T> structures;
    dataSet->getStructures(structures);
    const size_t numAtoms = dataSet->getNumAtomsWithRestriction();
    size_t idx = 0;
    for (const Structure<T> structure : structures.structures) {
        T base_inv[9];
        T* matrix = new T[numAtoms * numAtoms];
        std::fill_n(matrix, numAtoms * numAtoms, 0.0);
        DescriptorMatrix<T> descriptor;
        descriptor.numAtoms = structure.atoms.size();
        std::copy(structure.base, structure.base + 9, base_inv);
        lapack_int ipiv[3];
        lapack_int info = LapackGetrf<T>::call(
        /*LAPACK_ROW_MAJOR*/LAPACK_COL_MAJOR, 3, 3, base_inv, 3, ipiv);
        if (info < 0) {
            std::cerr << "Matrix factorization: " << "the " << -info
                    << "-th parameter of base matrix of structure " << idx
                    << " had an illegal value (skipping).";
            continue;
        } else if (info > 0) {
            std::cerr << "Matrix factorization: Base matrix of structure "
                    << idx << " is singular (skipping).";
            continue;
        }
        info = LapackGetri<T>::call(/*LAPACK_ROW_MAJOR*/LAPACK_COL_MAJOR, 3,
                base_inv, 3, ipiv);
        if (info < 0) {
            std::cerr << "Matrix invertation: " << "the " << -info
                    << "-th parameter of base matrix of structure " << idx
                    << " had an illegal value (skipping).";
            continue;
        } else if (info > 0) {
            std::cerr << "Matrix invertation: Base matrix of structure " << idx
                    << " is singular (skipping).";
            continue;
        }
        for (size_t i = 0; i < descriptor.numAtoms; i++) {
            const Atom<T> atom1 = structure.atoms[i];
            matrix[i * numAtoms + i] = 0.5 * pow(atom1.Z, 2.4);
            for (size_t j = i + 1; j < descriptor.numAtoms; j++) {
                const Atom<T> atom2 = structure.atoms[j];
                T r[3];
                MKLvSub<T>::call(3, atom1.r, atom2.r, r);
                T tmp[3] = { 0.0, 0.0, 0.0 };
                for (int k = 0; k < 3; k++) {
                    T x = 0.0;
                    for (int l = 0; l < 3; l++)
                        x += base_inv[3 * k + l] * r[l];
                    x = sin(x * M_PI);
                    x *= x;
                    for (int l = 0; l < 3; l++)
                        tmp[l] += structure.base[/*l * 3 + k*/k * 3 + l] * x;
                }
                T phi = 0.0;
                for (int k = 0; k < 3; k++)
                    phi += tmp[k] * tmp[k];
                phi = 1.0 / sqrt(phi);
                matrix[i * numAtoms + j] = atom1.Z * atom2.Z * phi;
                matrix[j * numAtoms + i] = matrix[i * numAtoms + j];
            }
        }
        this->sorting(matrix, numAtoms);
        descriptor.m = std::shared_ptr<T>(matrix);
        descriptor.name = structure.name;
        descriptors.push_back(descriptor);
        idx++;
    }
}

template class Descriptor<float> ;
template class Descriptor<double> ;
template class CoulombDescriptor<float> ;
template class SineDescriptor<double> ;
template class SineDescriptor<float> ;
