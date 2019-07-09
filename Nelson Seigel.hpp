//
//  Nelson Seigel.hpp
//  Grid Search For Nelson Seigel and yield Fit
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#ifndef Nelson_Seigel_hpp
#define Nelson_Seigel_hpp

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
using namespace std;
using vct_d=vector<double>;
using d_matrix=vector<vct_d>;

void datareader(string filename,map<double,double> &yield_data,map<double,double> &all_data);
void parser(stringstream& ss,double &time,double &yield);
void Grid_Search_NS(map<double,double> yield_data,map<double,double> all_data,vector<double> &parameters);
void NS_fit_yield_curve(vector<double> &parameters,map<double,double> &yield);
void dataoutput(string filename,map<double,double> &yield);

void BondReader(string file_name,map<string,vector<double>> &Bond_Info,string Valuation_date);
void parser_bond(stringstream & ss,string * str_data,int str_num);
void time_to_Maturity(vector<string> Maturity,vct_d &Accrual,
                      vct_d &time_to_maturity,string Valuate_date);
void parser_date(stringstream & ss,string * date_data);

void N_S_Bond_Price(map<string,vector<double>> Bond_Info,
                    map<string,vector<double>>& NewPrice,
                    vector<double> NS_parameters,string filename);

#endif /* Nelson_Seigel_hpp */
