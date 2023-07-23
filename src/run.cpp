#include "src/algs.hpp"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    auto A = init_diagonal_matrix_k1<double>(3,2);
    auto id = &A;
    auto dims = A.num_dimensions();
    long cols = A.shape()[1];
    auto shape = A.shape();
    cout<<"A="<<endl;
    print_t(cout, A);
    cout<<endl;
    auto I = inverse_diagonal_square_matrix(A[ boost::indices[range()<cols][range()] ]);
    id = &I;
    dims = I.num_dimensions();
    shape = I.shape();
    cout<<"I="<<endl;
    print_t(cout, I);
    cout<<endl;
    return 0;
}
