//
//  main.cpp
//  Vasicek_Model_Fitted_Yield
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#include <iostream>
#include<map>
#include<vector>
#include "Vasicek_Estimation1.hpp"
using namespace std;
int main(int argc, const char * argv[]) {
    //step 1: read in data
    string filename="//Users//daiyu//Desktop//TSA//Data Source//Yield_Curve.csv";
    map<double,double> yield_data;
    datareader(filename,yield_data);
    cout<<"***************************"<<endl;
    
    //Analytically result:
    vector<double> Vasicek_Parameters(4);
    Vasicek_estimation(yield_data,Vasicek_Parameters,4);
    
    //Fitted the Curve
    map<double,double> yield;
    Vasicek_fitted_Yield_Curve(Vasicek_Parameters,yield);
    
    //output file
    string filename1="//Users//daiyu//Desktop//TSA//Data Source//Vasicek_fitted_yield//Vasicek_Yield_Curve_1.csv";
    output_fitted_Yield_Curve(filename1,yield);
    
    // given term structure, find miss priced bond in the market
    //(1) read bond info
    string filename2="//Users//daiyu//Desktop//coupon bond YD -0411.csv";
    map<string,vector<double>> Bond_info;
    cout<<"***************************"<<endl;
    BondReader(filename2, Bond_info, "4/11/2019");
    
    //(2) price bond under Vasicek_term_structure
    map<string,vector<double>> NewPrice;
    
    string filename3="//Users//daiyu//Desktop//TSA//Data Source//Vasicek_fitted_yield//Vasicek_miss_Priced_Bond.csv";
    
    Vasicek_Bond_Price_1(Bond_info,NewPrice, Vasicek_Parameters,filename3);
    string filename4="//Users//daiyu//Desktop//TSA//Data Source//Vasicek_fitted_yield//Vasicek_Monte_carlo.csv";
    // Bond Price Evaluation Under  Simulation
    Vasicek_Simu1(Vasicek_Parameters,360,10,1.0/12,filename4);
    string filename5="//Users//daiyu//Desktop//TSA//Data Source//Vasicek_fitted_yield//Vasicek_Monte_carlo_2.csv";
    //Vasicek_Simu2(Vasicek_Parameters,360,10,1.0/12,filename5);
    return 0;
}
