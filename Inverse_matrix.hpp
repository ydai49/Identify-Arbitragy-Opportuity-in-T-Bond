//
//  Inverse_matrix.hpp
//  Term Structure From Par Bond
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#ifndef Inverse_matrix_hpp
#define Inverse_matrix_hpp

#include <stdio.h>
// REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>


void function_solver(boost::numeric::ublas::matrix<double> A,  boost::numeric::ublas::vector<double> B,boost::numeric::ublas::vector<double> &x)
{
    int n=(int)A.size1();
    boost::numeric::ublas::permutation_matrix<double> P(n);
    // fill matrix and rhs .. this is pretty straightforward
    boost::numeric::ublas::matrix<double> X(A);
    boost::numeric::ublas::vector<double> rhs(B);
    // fill matrix and rhs .. this is pretty straightforward

    // fill matrix and rhs .. this is pretty straightforward
    lu_factorize(X,P);
    // Now A and P contain the LU factorization of A
    x = rhs;
    lu_substitute(A,P,rhs);
}

#endif /* Inverse_matrix_hpp */
