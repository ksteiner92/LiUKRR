/*
 * types.h
 *
 *  Created on: 19.09.2016
 *      Author: Klaus Steiner
 */

#ifndef TYPES_H_
#define TYPES_H_

#include <mkl.h>
#include <matio.h>
#include <functional>
#include <map>
#include <utility>

typedef mat_t MatFileHandle;
typedef matvar_t MatVariable;

template<class F>
struct wrapper {
    static_assert(std::is_empty<F>(), "Lambdas must be empty");
    template<class ... Ts>
    decltype(auto) operator()(Ts&&... xs) const {
        return reinterpret_cast<const F&>(*this)(std::forward<Ts>(xs)...);
    }
};

template<class F>
constexpr wrapper<F> make_function(F*) {
    return {};
}

struct wrapper_factor {
    template<class F>
    constexpr wrapper<F> operator +=(F*) {
        return {};
    }
};

struct addr_add {
    template<class T>
    friend typename std::remove_reference<T>::type *operator+(addr_add, T &&t) {
        return &t;
    }
};

#define STATIC_LAMBDA wrapper_factor() += true ? nullptr : addr_add() + []

template<class T> struct Atom {
    T r[3];
    int Z;
};

template<class T> struct Structure {
    std::vector<Atom<T>> atoms;
    T base[9];
    T energy;
    std::string name;
};

template<class T> struct Structures {
    std::vector<Structure<T>> structures;
    size_t numAtoms;
};

struct SubSet {
    size_t *elements;
    size_t size;
    bool aligned;
};

template<class T> struct DescriptorMatrix {
    std::shared_ptr<T> m;
    int numAtoms;
    std::string name;
};

static std::map<std::string, int> PeriodicTable { { "H", 1 }, { "He", 2 }, {
        "Li", 3 }, { "Be", 4 }, { "B", 5 }, { "C", 6 }, { "N", 7 }, { "O", 8 },
        { "F", 9 }, { "Ne", 10 }, { "Na", 11 }, { "Mg", 12 }, { "Al", 13 }, {
                "Si", 14 }, { "P", 15 }, { "S", 16 }, { "Cl", 17 },
        { "Ar", 18 }, { "K", 19 }, { "Ca", 20 }, { "Sc", 21 }, { "Ti", 22 }, {
                "V", 23 }, { "Cr", 24 }, { "Mn", 25 }, { "Fe", 26 },
        { "Co", 27 }, { "Ni", 28 }, { "Cu", 29 }, { "Zn", 30 }, { "Ga", 31 }, {
                "Ge", 32 }, { "As", 33 }, { "Se", 34 }, { "Br", 35 },
        { "Kr", 36 }, { "Rb", 37 }, { "Sr", 38 }, { "Y", 39 }, { "Zr", 40 }, {
                "Nb", 41 }, { "Mo", 42 }, { "Tc", 43 }, { "Ru", 44 },
        { "Rh", 45 }, { "Pd", 46 }, { "Ag", 47 }, { "Cd", 48 }, { "In", 49 }, {
                "Sn", 50 }, { "Sb", 51 }, { "Te", 52 }, { "I", 53 },
        { "Xe", 54 }, { "Cs", 55 }, { "Ba", 56 }, { "La", 57 }, { "Ce", 58 }, {
                "Pr", 59 }, { "Nd", 60 }, { "Pm", 61 }, { "Sm", 62 },
        { "Eu", 63 }, { "Gd", 64 }, { "Tb", 65 }, { "Dy", 66 }, { "Ho", 67 }, {
                "Er", 68 }, { "Tm", 69 }, { "Yb", 70 }, { "Lu", 71 },
        { "Hf", 72 }, { "Ta", 73 }, { "W", 74 }, { "Re", 75 }, { "Os", 76 }, {
                "Ir", 77 }, { "Pt", 78 }, { "Au", 79 }, { "Hg", 80 },
        { "Tl", 81 }, { "Pb", 82 }, { "Bi", 83 }, { "Po", 84 }, { "At", 85 }, {
                "Rn", 86 }, { "Fr", 87 }, { "Ra", 88 }, { "Ac", 89 },
        { "Th", 90 }, { "Pa", 91 }, { "U", 92 }, { "Np", 93 }, { "Pu", 94 }, {
                "Am", 95 }, { "Cm", 96 }, { "Bk", 97 }, { "Cf", 98 },
        { "Es", 99 }, { "Fm", 100 }, { "Md", 101 }, { "No", 102 },
        { "Lr", 103 }, { "Rf", 104 }, { "Db", 105 }, { "Sg", 106 },
        { "Bh", 107 }, { "Hs", 108 }, { "Mt", 109 }, { "Ds", 110 },
        { "Rg", 111 }, { "Cn", 112 }, { "Uut", 113 }, { "Uuq", 114 }, { "Uup",
                115 }, { "Uuh", 116 }, { "Uus", 117 }, { "Uuo", 118 } };

/**
 * Template wrapper for cblas_?copy
 * Copies vector to another vector.
 */
template<class T> struct BlasCopy {
    static constexpr void (*call)(const MKL_INT n, const T *x,
            const MKL_INT incx, T *y, const MKL_INT incy) = NULL;
};
template<> struct BlasCopy<float> {
    static constexpr void (*call)(const MKL_INT n, const float *x,
            const MKL_INT incx, float *y, const MKL_INT incy) = &cblas_scopy;
};
template<> struct BlasCopy<double> {
    static constexpr void (*call)(const MKL_INT n, const double *x,
            const MKL_INT incx, double *y, const MKL_INT incy) = &cblas_dcopy;
};

/**
 * Template wrapper for cblas_?asum
 * Computes the sum of magnitudes of the vector elements.
 */
template<class T> struct BlasAsum {
    static constexpr T (*call)(const MKL_INT n, const T *x,
            const MKL_INT incx) = NULL;
};
template<> struct BlasAsum<float> {
    static constexpr float (*call)(const MKL_INT n, const float *x,
            const MKL_INT incx) = &cblas_sasum;
};
template<> struct BlasAsum<double> {
    static constexpr double (*call)(const MKL_INT n, const double *x,
            const MKL_INT incx) = &cblas_dasum;
};

template<class T> struct BlasNrm2 {
    static constexpr T (*call)(const MKL_INT n, const T *x,
            const MKL_INT incx) = NULL;
};
template<> struct BlasNrm2<float> {
    static constexpr float (*call)(const MKL_INT n, const float *x,
            const MKL_INT incx) = &cblas_snrm2;
};
template<> struct BlasNrm2<double> {
    static constexpr double (*call)(const MKL_INT n, const double *x,
            const MKL_INT incx) = &cblas_dnrm2;
};

/*
 * ( const char* uplo, const MKL_INT* n, const MKL_INT* nrhs, float* a,
 const MKL_INT* lda, MKL_INT* ipiv, float* b, const MKL_INT* ldb,
 float* work, const MKL_INT* lwork, MKL_INT* info )
 */
template<class T> struct Sysv {
    static constexpr void (*call)(const char* uplo, const MKL_INT* n,
            const MKL_INT* nrhs, T* a, const MKL_INT* lda, MKL_INT* ipiv, T* b,
            const MKL_INT* ldb, T* work, const MKL_INT* lwork,
            MKL_INT* info) = NULL;
};
template<> struct Sysv<float> {
    static constexpr void (*call)(const char* uplo, const MKL_INT* n,
            const MKL_INT* nrhs, float* a, const MKL_INT* lda, MKL_INT* ipiv,
            float* b, const MKL_INT* ldb, float* work, const MKL_INT* lwork,
            MKL_INT* info) = &ssysv;
};
template<> struct Sysv<double> {
    static constexpr void (*call)(const char* uplo, const MKL_INT* n,
            const MKL_INT* nrhs, double* a, const MKL_INT* lda, MKL_INT* ipiv,
            double* b, const MKL_INT* ldb, double* work, const MKL_INT* lwork,
            MKL_INT* info) = &dsysv;
};

/**
 * lapack_int LAPACKE_ssysv (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , float * a , lapack_int lda ,
 *  lapack_int * ipiv , float * b , lapack_int ldb );
 */
template<class T> struct LapackSysv {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, lapack_int nrhs, T * a, lapack_int lda,
    lapack_int * ipiv, T * b, lapack_int ldb) = NULL;
};
template<> struct LapackSysv<float> {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, lapack_int nrhs, float * a, lapack_int lda,
    lapack_int * ipiv, float * b, lapack_int ldb) = &LAPACKE_ssysv;
};
template<> struct LapackSysv<double> {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, lapack_int nrhs, double * a, lapack_int lda,
    lapack_int * ipiv, double * b, lapack_int ldb) = &LAPACKE_dsysv;
};

template<class T> struct MKLvSub {
    static constexpr void (*call)(const MKL_INT n, const T* a, const T* b,
            T *y) = NULL;
};
template<> struct MKLvSub<float> {
    static constexpr void (*call)(const MKL_INT n, const float* a,
            const float* b, float *y) = &vsSub;
};
template<> struct MKLvSub<double> {
    static constexpr void (*call)(const MKL_INT n, const double* a,
            const double* b, double *y) = &vdSub;
};

template<class T> struct MKLlanSy {
    static constexpr T (*call)(const char* norm, const char* uplo,
            const MKL_INT* n, const T* a, const MKL_INT* lda, T* work) = NULL;
};
template<> struct MKLlanSy<float> {
    static constexpr float (*call)(const char* norm, const char* uplo,
            const MKL_INT* n, const float* a, const MKL_INT* lda,
            float* work) = &slansy;
};
template<> struct MKLlanSy<double> {
    static constexpr double (*call)(const char* norm, const char* uplo,
            const MKL_INT* n, const double* a, const MKL_INT* lda,
            double* work) = &dlansy;
};

template<class T> struct LapackLansy {
    static constexpr T (*call)(int matrix_layout, char norm, char uplo,
    lapack_int n, const T * a, lapack_int lda) = NULL;
};
template<> struct LapackLansy<float> {
    static constexpr float (*call)(int matrix_layout, char norm, char uplo,
    lapack_int n, const float * a, lapack_int lda) = &LAPACKE_slansy;
};
template<> struct LapackLansy<double> {
    static constexpr double (*call)(int matrix_layout, char norm, char uplo,
    lapack_int n, const double * a, lapack_int lda) = &LAPACKE_dlansy;
};
/*
 * lapack_int LAPACKE_ssytri2 (int matrix_layout , char uplo , lapack_int n , float * a , lapack_int lda , const lapack_int * ipiv );
 */
template<class T> struct LapackSytri {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, T * a, lapack_int lda, const lapack_int * ipiv) = NULL;
};
template<> struct LapackSytri<float> {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, float * a, lapack_int lda,
            const lapack_int * ipiv) = &LAPACKE_ssytri;
};
template<> struct LapackSytri<double> {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, double * a, lapack_int lda,
            const lapack_int * ipiv) = &LAPACKE_dsytri;
};

/**
 * ssytrf
 * (int matrix_layout , char uplo , lapack_int n , float * a , lapack_int lda , lapack_int * ipiv )
 */
template<class T> struct LapackSsytrf {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, T * a, lapack_int lda, lapack_int * ipiv) = NULL;
};
template<> struct LapackSsytrf<float> {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, float * a, lapack_int lda,
    lapack_int * ipiv) = &LAPACKE_ssytrf;
};
template<> struct LapackSsytrf<double> {
    static constexpr lapack_int (*call)(int matrix_layout, char uplo,
    lapack_int n, double * a, lapack_int lda,
    lapack_int * ipiv) = &LAPACKE_dsytrf;
};

/*
 * void cblas_ssymv (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const float alpha, const float *a, const MKL_INT lda, const float *x, const MKL_INT incx, const float beta, float *y, const MKL_INT incy);
 *
 */
template<class T> struct BlasSsymv {
    static constexpr void (*call)(const CBLAS_LAYOUT Layout,
            const CBLAS_UPLO uplo, const MKL_INT n, const T alpha, const T *a,
            const MKL_INT lda, const T *x, const MKL_INT incx, const T beta,
            T *y, const MKL_INT incy) = NULL;
};
template<> struct BlasSsymv<float> {
    static constexpr void (*call)(const CBLAS_LAYOUT Layout,
            const CBLAS_UPLO uplo, const MKL_INT n, const float alpha,
            const float *a, const MKL_INT lda, const float *x,
            const MKL_INT incx, const float beta, float *y,
            const MKL_INT incy) = &cblas_ssymv;
};
template<> struct BlasSsymv<double> {
    static constexpr void (*call)(const CBLAS_LAYOUT Layout,
            const CBLAS_UPLO uplo, const MKL_INT n, const double alpha,
            const double *a, const MKL_INT lda, const double *x,
            const MKL_INT incx, const double beta, double *y,
            const MKL_INT incy) = &cblas_dsymv;
};

/**
 * lapack_int LAPACKE_sgetrf (int matrix_layout , lapack_int m , lapack_int n , float * a , lapack_int lda , lapack_int * ipiv );
 */
template<class T> struct LapackGetrf {
    static constexpr lapack_int (*call)(int matrix_layout, lapack_int m,
    lapack_int n, T * a, lapack_int lda, lapack_int * ipiv) = NULL;
};
template<> struct LapackGetrf<float> {
    static constexpr lapack_int (*call)(int matrix_layout, lapack_int m,
    lapack_int n, float * a, lapack_int lda,
    lapack_int * ipiv) = &LAPACKE_sgetrf;
};
template<> struct LapackGetrf<double> {
    static constexpr lapack_int (*call)(int matrix_layout, lapack_int m,
    lapack_int n, double * a, lapack_int lda,
    lapack_int * ipiv) = &LAPACKE_dgetrf;
};

template<class T> struct LapackGetri {
    static constexpr lapack_int (*call)(int matrix_layout, lapack_int n, T * a,
    lapack_int lda, const lapack_int * ipiv) = NULL;
};
template<> struct LapackGetri<float> {
    static constexpr lapack_int (*call)(int matrix_layout, lapack_int n,
            float * a,
            lapack_int lda, const lapack_int * ipiv) = &LAPACKE_sgetri;
};
template<> struct LapackGetri<double> {
    static constexpr lapack_int (*call)(int matrix_layout, lapack_int n,
            double * a, lapack_int lda,
            const lapack_int * ipiv) = &LAPACKE_dgetri;
};

template<class T> struct BlasGemv {
    static constexpr void (*call)(const CBLAS_LAYOUT Layout,
            const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n,
            const T alpha, const T *a, const MKL_INT lda, const T *x,
            const MKL_INT incx, const T beta, T *y, const MKL_INT incy) = NULL;
};
template<> struct BlasGemv<float> {
    static constexpr void (*call)(const CBLAS_LAYOUT Layout,
            const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n,
            const float alpha, const float *a, const MKL_INT lda,
            const float *x, const MKL_INT incx, const float beta, float *y,
            const MKL_INT incy) = &cblas_sgemv;
};
template<> struct BlasGemv<double> {
    static constexpr void (*call)(const CBLAS_LAYOUT Layout,
            const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n,
            const double alpha, const double *a, const MKL_INT lda,
            const double *x, const MKL_INT incx, const double beta, double *y,
            const MKL_INT incy) = &cblas_dgemv;
};

template<class T> struct ConvertToDoubleArray {
    static double* convert(T* in, const size_t n) {
        double* res = new double[n];
        for (size_t i = 0; i < n; i++)
            res[i] = in[i];
        delete[] in;
        return res;
    }
};
template<> struct ConvertToDoubleArray<double> {
    static double* convert(double* in, const size_t n) {
        return in;
    }
};

#endif /* TYPES_H_ */
