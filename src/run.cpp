#include "src/algs.hpp"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    auto A = init_triangular_matrix_k1<int>(4,2);
    auto dims = A.num_dimensions();
    long cols = A.shape()[1];
    auto shape = A.shape();
    cout<<"A="<<endl;
    print_t(cout, A);
    cout<<endl;
    auto I = inverse_triangular_square_matrix(A[ boost::indices[range()<cols][range()] ]);
    dims = I.num_dimensions();
    shape = I.shape();
    cout<<"I="<<endl;
    print_t(cout, I);
    cout<<endl;
    return 0;
}
