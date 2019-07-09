//
//  Par Bond.hpp
//  Term Structure From Par Bond
//
//  Created by 代雨 on 2019/5/1.
//  Copyright © 2019年 代雨. All rights reserved.
//

#ifndef Par_Bond_hpp
#define Par_Bond_hpp

#include <stdio.h>
#include<vector>
#include<map>
#include<iostream>
#include<string>
#include<boost/date_time/gregorian/gregorian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace std;
using vct_d=vector<double>;
using d_matrix=vector<vct_d>;

d_matrix cleaned_data(string file_name,string Valuation_date);
void BondReader(string file_name,map<string,vector<double>> &Bond_Info,string Valuation_date);

void parser(stringstream & ss,string * str_data,int str_num);
void Yield_Curve_Output(string filename,map<double,double> Linear_Yield_Curve);
void time_to_Maturity(vector<string> Maturity,vct_d &Accrual,
                       vct_d &time_to_maturity,string Valuate_date);
void parser_date(stringstream & ss,string * date_data);
double N_S_Method_Yield(double coupon,double quoted_price,double accrual,double time_to_maturity);
double quoted_price_1st_dev(double coupon,double accrual,double
                            time_to_maturity,double y);
double Quoted_Price(double coupon,double price,double accrual,double
                    time_to_maturity,double y);
void sopt_rate_from_YTM();
void linear_intropolate_YTM(map<double,double> time_yield,
                            map<double,double> &time_yield_intropolate,
                            vector<double> TTM);
void Cubic_Boostrapping(map<double,double> time_yield,
                        map<double,double> &cubic_bootstrap,
                        vector<double> TTM);
void fun_solver(boost::numeric::ublas::matrix<double> &A,  boost::numeric::ublas::vector<double> &B,boost::numeric::ublas::vector<double> &x);

//Pricing Bond
void BondReader(string file_name,map<string,vector<double>> &Bond_Info,string Valuation_date);

#endif /* Par_Bond_hpp */
