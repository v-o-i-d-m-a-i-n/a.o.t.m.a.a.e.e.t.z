#pragma once

#include <iostream>
#include <assert.h>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/next_prior.hpp>

using range = boost::multi_array_types::index_range;

template <typename T, std::size_t dims>
using tensor = boost::multi_array<T,dims>;

template <typename T, std::size_t dims>
using tensor_view = boost::detail::multi_array::multi_array_view<T,dims>;

template <typename T>
using matrix = tensor<T,2>;

template <typename T>
using matrix_view = tensor_view<T,2>;

auto& extents = boost::extents;

template <typename T, std::size_t num>
using array = boost::array<T,num>;

template <std::size_t dims>
using index = array<std::size_t,dims>;

template <typename T, std::size_t dims>
void print_t(std::ostream& os, const tensor<T,dims>& A, char method = 1);

template <typename T, std::size_t dims>
tensor<T, dims> set_value(tensor<T, dims> A, T value);

template <typename T>
matrix<T> init_triangular_matrix_k1(const std::size_t& num_rows, const std::size_t& k);

template <typename Tin, typename Tout>
matrix<Tout> inverse_triangular_square_matrix(const matrix<Tin>& A);

template <typename Tin, typename Tout>
matrix<Tout> inverse_triangular_square_matrix(const matrix_view<Tin>& A);

template <typename T>
matrix<double> inverse_triangular_square_matrix_real(const matrix<T>& A);

template <typename T>
const char triangular_matrix_assertion(const matrix<T>& mat);

template <typename T, std::size_t dims>
void print_t(std::ostream& os, const tensor<T,1>& A, char method = 1)
{
    if (method==1) {
        using boost::next;
        typename tensor<T,dims>::const_iterator i;
        
        os << "[";
        for (i = A.begin(); i != A.end(); ++i) {
            //print(os, *i); //
            os << *i;
            if (next(i) != A.end()) {
                os << ',';
            }
        }
        os << "]";
    } else {
        auto ndims = A.num_dimensions();
        auto shape = A.shape();
        index<dims> idx = {};
        std::size_t i=0, _i;
        if (ndims>0) {
            for (_i=0; _i<ndims; ++_i) {
                os << "[";
            }
            _i -= 1;
            while(true){
                if (idx[_i]<shape[_i]) {
                    if (i>0) {
                        os << std::endl;
                        for (std::size_t k=0; k<ndims; ++k) {
                            if (k<ndims-i) {
                                os << " ";
                            } else {
                                os << "[";
                            }
                        }
                    }
                    os << A(idx) << " ";
                    i=0;
                } else {
                    os << "]";
                    idx[_i]=0;
                    ++i;
                }
                if (i<ndims) {
                    _i = ndims-1-i;
                    idx[_i]++;
                } else {
                    break;
                }
            }
        }
    }
}

template <typename T, std::size_t dims>
void print_t(std::ostream& os, const tensor<T,dims>& A, char method)
{
    if (method==1) {
        using boost::next;
        typename tensor<T,dims>::const_iterator i;
        
        os << "[" ;
        for (i = A.begin(); i != A.end(); ++i) {
            print_t<T,dims-1>(os, *i, 1);
            if (next(i) != A.end()){
                os << ',' << std::endl;
            }
        }
        os << "]";
    } else {
        auto ndims = A.num_dimensions();
        auto shape = A.shape();
        index<dims> idx = {};
        std::size_t i=0, _i;
        if (ndims>0) {
            for (_i=0; _i<ndims; ++_i) {
                os << "[";
            }
            _i -= 1;
            while(true){
                if (idx[_i]<shape[_i]) {
                    if (i>0) {
                        os << std::endl;
                        for (std::size_t k=0; k<ndims; ++k) {
                            if (k<ndims-i) {
                                os << " ";
                            } else {
                                os << "[";
                            }
                        }
                    }
                    os << A(idx) << " ";
                    i=0;
                } else {
                    os << "]";
                    idx[_i]=0;
                    ++i;
                }
                if (i<ndims) {
                    _i = ndims-1-i;
                    idx[_i]++;
                } else {
                    break;
                }
            }
        }
    }
}

template <typename T, std::size_t dims>
tensor<T, dims> set_value(tensor<T, dims> A, T value){
    auto ndims = A.num_dimensions();
    auto shape = A.shape();
    index<dims> idx = {};
    A(idx) = value;
    auto ele = A(idx); //dbg
    std::size_t ii = 0;
    while(ii<ndims){
        if (idx[ii]<shape[ii]-1) {
            idx[ii]++;
            ele = A(idx);//dbg
            A(idx) = value;
            ele = A(idx);//dbg
            ii=0;
        } else {
            idx[ii]=0;
            ii++;
        }
    }
    return A;
}

template <typename T>
matrix<T> init_triangular_matrix_k1(const std::size_t& num_rows, const std::size_t& k)
{
    std::size_t nr = num_rows, nc = num_rows-k+1;
    matrix<T> A(extents[nr][nc]);
    for (std::size_t j=0;j<nc;++j) {
        for (std::size_t i=j;i<j+k;++i) {
            A[i][j]=T(-1);
        }
    }
    return A;
}

matrix<double> inverse_triangular_square_matrix(const matrix<double>& A){
    return inverse_triangular_square_matrix_real(A);
}

matrix<double> inverse_triangular_square_matrix(const matrix_view<double>& A)
{
    // matrix<double> Inv = inverse_triangular_square_matrix(matrix<double>(A)); 
    return inverse_triangular_square_matrix(matrix<double>(A));
}

matrix<double> inverse_triangular_square_matrix(const matrix<int>& A){
    return inverse_triangular_square_matrix_real(A);
}

matrix<double> inverse_triangular_square_matrix(const matrix_view<int>& A)
{
    return inverse_triangular_square_matrix(matrix<int>(A));
}

template <typename T>
matrix<double> inverse_triangular_square_matrix_real(const matrix<T>& A)
{
    auto shape = A.shape();
    long nr = shape[0]; long nc = shape[1];
    assert(nr==nc);
    auto LU_detect = triangular_matrix_assertion(A);
    matrix<double> LU(extents[nr][nc]);
    matrix<double> Inv(extents[nr][nc]);
    
    // Step 1, a(ik) = a(ik)/a(kk) k ∈ [0,i-1]
    for (std::size_t i=0; i<nr; ++i) {
        assert(A[i][i]!=T(0));
        LU[i][i] = A[i][i];
        for (std::size_t k=0; k<i; ++k) {
            LU[i][k] = A[i][k]; LU[i][k] /= LU[k][k];
        }
    }
    
    if (LU_detect=='U'){
    //    LU = permute(LU);
    }
    
    for (std::size_t i=0; i<nr; ++i) {
        // Step 2, b(ii) = 1
        Inv[i][i]=1.0;
        // Step 3, b(ij) = -SUM(k=j to i-1)(a(ik)*b(kj)) (j<i)
        for (std::size_t j=0; j<i; ++j) {
            for (std::size_t k=j; k<=i-1; ++k) {
                Inv[i][j] -= LU[i][k] * Inv[k][j];
            } 
        }
    }

    // Step 4, b(ij) = b(ij)/a(ii)
    for (std::size_t i=0; i<nr; ++i) {
        for (std::size_t j=0; j<=i; ++j) {
            Inv[i][j] /= LU[i][i];
        }
    }

    if (LU_detect=='U'){
    //    Inv = permute(Inv);
    }
    
    return Inv;
}

template <typename T>
const char triangular_matrix_assertion(const matrix<T>& mat){
    auto shape = mat.shape();
    long nr = shape[0]; 
    long nc = shape[1];
    char LU_detect = 'D';
    
    for (std::size_t i=0; i<nr; ++i) {
        for (std::size_t j=0; j<nc; ++j) {
            if (mat[i][j]!=T(0)) {
                if (i==j) {
                    LU_detect = LU_detect;
                } else {
                    if (i>j) {
                        assert(LU_detect!='U' && "The matrix is not triangle.");
                        LU_detect = 'L';
                    } else {
                        assert(LU_detect!='L' && "The matrix is not triangle.");
                        LU_detect = 'U';
                    }
                }
            }
        }
    };
    
    return LU_detect;
}
