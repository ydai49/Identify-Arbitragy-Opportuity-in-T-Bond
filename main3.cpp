//
//  main.cpp
//  Grid Search For Nelson Seigel and yield Fit
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#include <iostream>
#include<map>
#include<vector>

#include"Nelson Seigel.hpp"
int main(int argc, const char * argv[]) {
    
    //step 1: choosing data point to do estimation
    string filename="//Users//daiyu//Desktop//TSA//Data Source//Yield_Curve.csv";
    map<double,double> yield_data;
    map<double,double> all_data;
    datareader(filename,yield_data,all_data);
    cout<<"***************************"<<endl;
   
    // Step 2: Apply grid search method using OLS to find parameters for NS_model
    vector<double> NS_parameters(4);
    Grid_Search_NS(yield_data,all_data,NS_parameters);
    cout<<"***************************"<<endl;

    // Step 3: Fit the yield with NS model
    map<double,double> yield;
    NS_fit_yield_curve(NS_parameters,yield);

    // Step 4: output yield to csv
    string filename1="//Users//daiyu//Desktop//TSA//Data Source//NS_grid_fit//NS_Yield_Curve_fit.csv";
    dataoutput(filename1,yield);

    //Step 5: given term structure, find miss priced bond in the market
    //(1) read bond indo
    string filename3="//Users//daiyu//Desktop//coupon bond YD -0411.csv";
    map<string,vector<double>> Bond_info;
    BondReader(filename3, Bond_info, "4/11/2019");
    //(2) price bond under NS_term_structure
    map<string,vector<double>> NewPrice;
    string filename4="//Users//daiyu//Desktop//TSA//Data Source//NS_grid_fit//NS_miss_Priced_Bond.csv";
    N_S_Bond_Price(Bond_info,NewPrice, NS_parameters,filename4);
    return 0;
}
