//
//  Vasicek_Estimation1.hpp
//  Vasicek_Model_Fitted_Yield
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#ifndef Vasicek_Estimation1_hpp
#define Vasicek_Estimation1_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>
#include<string>
#include<vector>
#include<map>
#include<sstream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<boost/date_time/gregorian/gregorian.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include<Python/Python.h>
using namespace std;
using vct_d=vector<double>;
using d_matrix=vector<vct_d>;

void datareader(string filename,map<double,double> &yield_data);
void parser(stringstream& ss,double &time,double &yield);
void Vasicek_estimation(map<double,double> yield_data,vector<double>& result,int num_of_result);
void Vasicek_fitted_Yield_Curve(vector<double> result,map<double,double>& yield);
void output_fitted_Yield_Curve(string filename,map<double,double> yield);

void BondReader(string file_name,map<string,vector<double>> &Bond_Info,string Valuation_date);
void parser_bond(stringstream & ss,string * str_data,int str_num);
void time_to_Maturity(vector<string> Maturity,vct_d &Accrual,
                      vct_d &time_to_maturity,string Valuate_date);
void parser_date(stringstream & ss,string * date_data);
void Vasicek_Bond_Price_1(map<string,vector<double>> Bond_Info,
                     map<string,vector<double>> &NewPrice,
                     vector<double> Vasicek_parameters,string filename);
//this is a function that use Box-Muller Method to generate random number
double BoxMullerMethod();
void Vasicek_Simu1(vector<double> Vasicek_parameters,int num,int path,double dt,string filename);
//void Vasicek_Simu2(vector<double> Vasicek_parameters,int num,int path,double dt,string filename);
#endif /* Vasicek_Estimation1_hpp */
