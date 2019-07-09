//
//  main.cpp
//  Term Structure From Par Bond
//
//  Created by 代雨 on 2019/5/1.
//  Copyright © 2019年 代雨. All rights reserved.
//

#include <iostream>
#include<vector>
#include<algorithm>
#include<fstream>
#include<iterator>
#include<map>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "Par Bond.hpp"
using namespace std;

using vct_d=vector<double>;
using d_matrix=vector<vct_d>;
int main(int argc, const char * argv[]) {
    // Calculate Yield From Par Bond and retrieve spot rate from par bond yield
    string filename="//Users//daiyu//Desktop//TSA//Data Source//coupon bond YD -0411.csv";
    d_matrix my_bond_info=cleaned_data(filename,"4/11/2019");
    int row=(int) my_bond_info[0].size();
    //calculated YTM from par bond
    vector<double> YTM;
    vector<double> TTM; //time to maturity
    map<double, double> Yield_Curve;
    map<double,double> Linear_Yield_Curve;
    
    for(int i=0;i<row;i++){
        YTM.push_back(N_S_Method_Yield(my_bond_info[0][i], my_bond_info[1][i], my_bond_info[2][i], my_bond_info[3][i]));
        TTM.push_back(my_bond_info[2][i]);
        Yield_Curve[TTM[i]]=YTM[i];
    }
    
    string yield_curve="//Users//daiyu//Desktop//TSA//Data Source//Par_Bond//par_bond_yield.csv";
    Yield_Curve_Output(yield_curve,Yield_Curve);
    
    //Cook Term structure from Yield Curve
    //Method 1 Linear Boostrapping
    linear_intropolate_YTM(Yield_Curve,Linear_Yield_Curve,TTM);
    string Linear_Boostrapping="//Users//daiyu//Desktop//TSA//Data Source//Par_Bond//Linear_Yield_Curve.csv";
    Yield_Curve_Output(Linear_Boostrapping,Linear_Yield_Curve);
    
    //Method 2 Cubic Boostrapping
    map<double,double> cubic_bootstrap;
    Cubic_Boostrapping(Yield_Curve,cubic_bootstrap,TTM);
    string Cubic_Boostrapping="//Users//daiyu//Desktop//TSA//Data Source//Par_Bond//Cubic_Yield_Curve.csv";
    Yield_Curve_Output(Cubic_Boostrapping,cubic_bootstrap);
    
    return 0;
}
