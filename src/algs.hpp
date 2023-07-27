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
void print_t(std::ostream& os, const tensor<T,dims>& A);

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
void print_t(std::ostream& os, const tensor<T,1>& A)
{
    using boost::next;
    typename tensor<T,dims>::const_iterator i;
    
    os << "[";
    for (i = A.begin(); i != A.end(); ++i) {
        //print(os, *i); //
        os << *i;
        if (next(i) != A.end()){
            os << ',';
        }
    }
    os << "]";
}

template <typename T, std::size_t dims>
void print_t(std::ostream& os, const tensor<T,dims>& A)
{
    using boost::next;
    typename tensor<T,dims>::const_iterator i;
    
    os << "[" ;
    for (i = A.begin(); i != A.end(); ++i) {
        print_t<T,dims-1>(os, *i);
        if (next(i) != A.end()){
            os << ',' << std::endl;
        }
    }
    os << "]";
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
        if (idx[ii]<shape[ii]-1){
            idx[ii]++;
            ele = A(idx);//dbg
            A(idx) = value;
            ele = A(idx);//dbg
            ii=0;
        }else{
            idx[ii]=0;
            ii++;
        }
    }
    return A;
}

template <typename T>
matrix<T> init_triangular_matrix_k1(const std::size_t& num_rows, const std::size_t& k)
{
	std::size_t nr = num_rows, nc = num_rows-k;
	matrix<T> A(extents[nr][nc]);
	for(std::size_t i=0;i<nr;++i){
		for(std::size_t j=0;j<nc&&j<=i;++j){
			A[i][j]=1;
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

template <typename T>
matrix<double> inverse_triangular_square_matrix_real(const matrix<T>& A)
{
    auto shape = A.shape();
    long nr = shape[0]; long nc = shape[1];
    assert(nr==nc);
    auto LU_detect = triangular_matrix_assertion(A);
    matrix<double> LU(A);
    matrix<double> Inv(extents[nr][nc]);

    if (LU_detect=='U'){
    //    LU = permute(LU);
    }

    // Step 1, a(ik) = a(ik)/a(kk)/a(ii) k ∈ [0,i-1]
    for (long i=0; i<nr; ++i) {
        for (long k=0; k<i; ++k) {
            assert(LU[k][k]*LU[i][i]!=0.0);
            LU[i][k] /= (LU[k][k]*LU[i][i]);
        }
    }
    for (long i=0; i<nr; ++i) {
        // Step 2, b(ii) = 1
        Inv[i][i]=1.0;
        // Step 3, b(ij) = -SUM(k=j to i-1)(a(ik)*b(kj)) (j<i)
        for (long j=0; j<i; ++j) {
            for (long k=j; k<=i-1; ++k) {
                Inv[i][j] -= LU[i][k] * Inv[k][j];
            }
        }
        // Step 4, b(ii) = b(ii)/a(ii)
        Inv[i][i]=1.0/A[i][i];
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
    
    for (long i=0; i<nr; ++i) {
        for (long j=0; j<nc; ++j) {
            if (mat[i][j]!=T(0)) {
                if (i==j) {
                    LU_detect = LU_detect;
                }else {
                    if (i>j) {
                        assert(LU_detect!='U' && "The matrix is not triangle.");
                        LU_detect = 'L';
                    }
                    else {
                        assert(LU_detect!='L' && "The matrix is not triangle.");
                        LU_detect = 'U';
                    }
                }
            }
        }
    };
    
    return LU_detect;
}
