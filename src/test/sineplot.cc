/*
 * sineplot.cc
 *
 *  Created on: 12.11.2016
 *      Author: Klaus Steiner
 */

#include <iostream>
#include <fstream>
#include <random>

#include "machine.h"
#include "dataset.h"
#include "descriptor.h"
#include "progargparser.h"

const double dot(const double *a, const double *b) {
    double dot = 0.0;
    for (int i = 0; i < 3; i++)
        dot += a[i] * b[i];
    return dot;
}

int main(int argc, char** argv) {
    VASPDataSet *dataSet = new VASPDataSet(
            "../../data/mp-basicStructures_Symmetrized/", true);
    Structures<double> structures;
    std::uniform_int_distribution<size_t> dist(0,
            dataSet->getTotalNumCompounds() - 1);
    std::random_device rd; // uses RDRND or /dev/urandom
    //const size_t idx = dist(rd);
    const size_t idx = 100;
    dataSet->setSubSet(idx, idx + 1);
    dataSet->getStructures(structures);
    const size_t numAtoms = structures.numAtoms;
    const Structure<double> structure = structures.structures[0];
    double base_inv[9];
    std::copy(structure.base, structure.base + 9, base_inv);
    lapack_int ipiv[3];
    lapack_int info = LapackGetrf<double>::call(
    LAPACK_ROW_MAJOR/*LAPACK_COL_MAJOR*/, 3, 3, base_inv, 3, ipiv);
    if (info < 0) {
        std::cerr << "Matrix factorization: " << "the " << -info
                << "-th parameter of base matrix of structure " << idx
                << " had an illegal value (skipping).";
        return 1;
    } else if (info > 0) {
        std::cerr << "Matrix factorization: Base matrix of structure " << idx
                << " is singular (skipping).";
        return 1;
    }
    info = LapackGetri<double>::call(LAPACK_ROW_MAJOR/*LAPACK_COL_MAJOR*/, 3,
            base_inv, 3, ipiv);
    if (info < 0) {
        std::cerr << "Matrix invertation: " << "the " << -info
                << "-th parameter of base matrix of structure " << idx
                << " had an illegal value (skipping).";
        return 1;
    } else if (info > 0) {
        std::cerr << "Matrix invertation: Base matrix of structure " << idx
                << " is singular (skipping).";
        return 1;
    }
    double r1[3] = { 0.0, 0.0, 0.0 };
    double r2[3] = { 0.0, 0.0, 0.0 };
    double x[3];
    double y[3];
    double incx[3];
    double incy[3];
    std::copy(structure.base, structure.base + 3, x);
    std::copy(structure.base + 3, structure.base + 6, y);
    double tmp[3];
    for (int i = 0; i < 3; i++)
        tmp[i] = y[i] - dot(x, y) / dot(x, x) * x[i];
    std::copy(tmp, tmp + 3, y);
    for (int i = 0; i < 3; i++) {
        incx[i] = x[i] / 50;
        incy[i] = y[i] / 50;
    }
    std::ofstream fs;
    fs.open("/home/klaus/dev/LiU/test/sineplot.dat", std::ios::out);
    fs << "# x" << std::setw(30) << "y" << std::setw(30) << "val" << std::endl;
    for (size_t i = 0; i <= 200; i++) {
        for (int k = 0; k < 3; k++)
            r1[k] += incx[k];
        std::fill_n(r2, 3, 0.0);
        for (size_t j = 0; j <= 200; j++) {
            for (int k = 0; k < 3; k++)
                r2[k] += incy[k];
            double r[3];
            MKLvSub<double>::call(3, r1, r2, r);
            double tmp[3] = { 0.0, 0.0, 0.0 };
            for (int k = 0; k < 3; k++) {
                double x = 0.0;
                for (int l = 0; l < 3; l++)
                    x += base_inv[3 * k + l] * r[l];
                x = sin(x * M_PI);
                x *= x;
                for (int l = 0; l < 3; l++)
                    tmp[l] += structure.base[l * 3 + k/*k * 3 + l*/] * x;
            }
            double phi = 0.0;
            for (int k = 0; k < 3; k++)
                phi += tmp[k] * tmp[k];
            phi = phi > 1e-6 ? 1.0 / sqrt(phi) : 0.1;
            //phi = 1.0 / sqrt(phi);
            fs << i << std::setw(30) << j << std::setw(30) << phi << std::endl;
        }
        fs << std::endl;
    }
    fs.close();
    delete dataSet;
    return 0;
}

