//
//  Par Bond.cpp
//  Term Structure From Par Bond
//
//  Created by 代雨 on 2019/5/1.
//  Copyright © 2019年 代雨. All rights reserved.
//

#include "Par Bond.hpp"
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<algorithm>
#include<string>
#include<map>
#include<math.h>
#include<sstream>
#include<boost/date_time/gregorian/gregorian.hpp>
#include"Inverse_matrix.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
using namespace std;
using vct_d=vector<double>;
using d_matrix=vector<vct_d>;

d_matrix cleaned_data(string file_name,string Valuation_date){
    // read in data
    ifstream infile;
    
    infile.open(file_name);
    if(infile.fail()){
        cerr<<"Can't open file"<<endl;
        exit(2);
    }
    vector<double> Coupon;
    vector<string> Maturity;
    vector<double> Ask_price;
    
    d_matrix Bond_Info;
    //read in header
    string readline;
    getline(infile,readline);

    //read in body,store all the information in bond_info as string type data
    while(getline(infile,readline)){
        stringstream ss(readline),new_ss;
        int Strnum=11;
        string str_data[Strnum];
        double coupon,price;
        
        //parser data //consider wrong data ,how will you handle it
        parser(ss,str_data,Strnum);
        new_ss<<str_data[2]<<'\t'<<str_data[9];
        new_ss>>coupon>>price;
        
        //only want par bond
        if(price<101 && price>99){
            Maturity.push_back(str_data[3]);
            Coupon.push_back(coupon);
            Ask_price.push_back(price);
        }
    }
    infile.close();
    
    int par_bond_num=int(Ask_price.size());
    vct_d Accrual(par_bond_num);
    vct_d time_to_maturity(par_bond_num);
    cout<<"Valuation_date "<<Valuation_date<<endl;

    time_to_Maturity(Maturity,Accrual,time_to_maturity,Valuation_date);
    Bond_Info.push_back(Coupon);
    Bond_Info.push_back(Ask_price);
    Bond_Info.push_back(time_to_maturity);
    Bond_Info.push_back(Accrual);
    return Bond_Info;
}

//*****************************************************************************//

void parser(stringstream & ss,string * str_data,int str_num)
{
    for(int i=0;i<str_num;i++)
        getline(ss,str_data[i],',');
}
//*****************************************************************************//
void parser_date(stringstream & ss,string * date_data)
{
    //4/11/2019 month/day/year
    for(int i=0;i<3;i++)
        getline(ss,date_data[i],'/');
}
//*****************************************************************************//
void time_to_Maturity(vector<string> Maturity,vct_d &Accrual,
                       vct_d &time_to_maturity,string Valuation_date)
{
    int Val_date[3];
    string those_days[3];
    stringstream ss(Valuation_date),new_ss;
    
    parser_date(ss,those_days);
    new_ss<<those_days[0]<<'\t'<<those_days[1]<<'\t'<<those_days[2];
    new_ss>>Val_date[0]>>Val_date[1]>>Val_date[2];
    boost::gregorian::date MY_Valuation_date{Val_date[2],Val_date[0],Val_date[1]};
    
    int count=0;
    //get time to maturity
    int Diff_Day[int(time_to_maturity.size())];
    for(auto &m_date :Maturity){
        stringstream ssr(m_date),new_ssr;
        parser_date(ssr,those_days);
        new_ssr<<those_days[0]<<'\t'<<those_days[1]<<'\t'<<those_days[2];
        new_ssr>>Val_date[0]>>Val_date[1]>>Val_date[2];
        boost::gregorian::date Bond_Maturity{Val_date[2],Val_date[0],Val_date[1]};
        boost::gregorian::days diff_days=Bond_Maturity-MY_Valuation_date;
        time_to_maturity[count]=diff_days.days()*1.0/365;
        Diff_Day[count]=int(diff_days.days());
        count++;
    }
    
    //get Accrual
    int basis=182;
    // Maturity-Last_coupon mod basis=integer
    // Maturity-now=(integer-1)+decimal_part
    //decimal_part=(basis-accural)
    for(int i=0;i<time_to_maturity.size();i++){
        int integer_part=floor(Diff_Day[i]/basis);
        int decimal_part=Diff_Day[i]-integer_part*basis;
        Accrual[i]=basis-decimal_part;
        //cout<<Accrual[i]<<" days\n";
    }
}
//*****************************************************************************//
// one important feature of par bond is that its YTM=coupon rate
double N_S_Method_Yield(double coupon,double quoted_price,double accrual,double time_to_maturity)
{
    double guess=coupon/100;
    double h =Quoted_Price(coupon,quoted_price,accrual,time_to_maturity,guess)/
        quoted_price_1st_dev(coupon,accrual,time_to_maturity,guess);
    double precision=0.00001;
    while(abs(h)>precision){
        guess=guess-h;
        h= Quoted_Price(coupon,quoted_price,accrual,time_to_maturity,guess)/
        quoted_price_1st_dev(coupon,accrual,time_to_maturity,guess);
    }
    //cout<<"the value of root is "<<guess<<endl;
    return guess;
}

double Quoted_Price(double coupon,double price,double accrual,double
                    time_to_maturity,double y)
{
    int basis=182;
    y=y/2;//semianual
    coupon=coupon/2;//semianual
    double a=(basis-accrual)/basis-1;
    double b=coupon*accrual*1.0/basis;
    int N=floor(time_to_maturity/0.5);
    double quoted_price=0.;
    double Coupon_part=(coupon)/(y*pow(1+y, a))*(1-1/pow(1+y, N));
    double Face_value=100/(pow(1+y, N+a));
    quoted_price=Coupon_part+Face_value-b-price;
    return quoted_price;
}

double quoted_price_1st_dev(double coupon,double accrual,double
                    time_to_maturity,double y)
{
    int basis=182;
    y=y/2;//semianual
    coupon=coupon/2;//semianual
    double a=(basis-accrual)/basis-1;
    double b=coupon*accrual*1.0/basis;
    int N=floor(time_to_maturity/0.5);
    double quoted_price=0.;
    
    double part1=(-1)*(1-1/pow(1+y, N))*coupon*(1/(y*y*pow(1+y,a))+a/
                                                   (y*pow(1+y,a+1)));
    double part2=coupon/(y*pow(1+y, a))*1/pow(1+y, N+1);
    double part3=100*(N+a)/pow(1+y, N+a+1);
    quoted_price=part1+part2-part3;
    return quoted_price;
}

void linear_intropolate_YTM(map<double,double> time_yield,
                                map<double,double> &time_yield_intropolate,
                                vector<double> TTM)
{
    double day_shift=1.0/365;
    vector<double> new_TTM(TTM);
    sort(new_TTM.begin(), new_TTM.end());
    for(int i=0;i<TTM.size()-1;i++){
        
        time_yield_intropolate[new_TTM[i]]=time_yield[new_TTM[i]];
        
        double new_time=new_TTM[i];
        
        int n=floor((new_TTM[i+1]-new_TTM[i])/day_shift);
        
        double d_yield=(time_yield[new_TTM[i+1]]-time_yield[new_TTM[i]])/n;
        
        for(int j=1;j<n;j++){
            time_yield_intropolate[new_TTM[i]+j*day_shift]=
                d_yield*j+time_yield[new_TTM[i]];
        }
    }

}
void Yield_Curve_Output(string filename,map<double,double> Linear_Yield_Curve)
{
    ofstream outfile;
    outfile.open(filename);
    outfile<<"time to maturity"<<'\t'<<"yield"<<endl;
    for(auto & e: Linear_Yield_Curve)
        outfile<<e.first<<'\t'<<e.second<<endl;
}
//Cubic Boostrapping
//yt=a*t^3+b*t^2+c*t+d
//we have yt and t==>solve for a,b,c,d
//In the matrix form : parameters=inv(t(X)*X)*t(X)*Y
void Cubic_Boostrapping(map<double,double> time_yield,
                        map<double,double> &cubic_bootstrap,
                        vector<double> TTM)
{
    vector<double> new_TTM(TTM);
    sort(new_TTM.begin(), new_TTM.end());
    //time and yield
    int n=(int)new_TTM.size();
    boost::numeric::ublas::matrix<double> X(n,4);
    boost::numeric::ublas::vector<double> Y(n);
    boost::numeric::ublas::vector<double> new_Y(4);
    for(int i=0;i<n;++i){
        Y(i)=time_yield[new_TTM[i]];
        for(int j=0;j<4;++j)
            X(i,j) = pow(new_TTM[i],3-j);
    }
    boost::numeric::ublas::matrix<double> A;
    boost::numeric::ublas::vector<double> parameters;
    A=prod(trans(X),X);
    new_Y=prod(trans(X),Y);
    boost::numeric::ublas::permutation_matrix<double> P(4);
    lu_factorize(A,P);
    // Now A and P contain the LU factorization of A
    parameters = new_Y;
    lu_substitute(A,P,parameters);
    cout<<"parameters: "<<parameters<<endl;//
    
    //Construct Yield Curve
    //yt=a*t^3+b*t^2+c*t+d
    double day_shift=1.0/365;
    for(int i=0;i<TTM.size()-1;i++){
        int n=floor((new_TTM[i+1]-new_TTM[i])/day_shift);
        for(int j=1;j<n;j++){
            double t=new_TTM[i]+j*day_shift;
            cubic_bootstrap[t]=parameters[0]*t*t*t+parameters[1]*t*t+
            parameters[2]*t+parameters[3];
        }
    }
}

//pricing Bond
void BondReader(string file_name,map<string,vector<double>> &Bond_Info,string Valuation_date){
    // read in data
    ifstream infile;
    
    infile.open(file_name);
    if(infile.fail()){
        cerr<<"Can't open file"<<endl;
        exit(2);
    }
    
    vector<string> Identifier;
    vector<string> Maturity;
    vector<double> T_t_M;//time to maturity
    vector<double> Ask_price;
    vector<double> Coupon;
    
    
    //read in header
    string readline;
    getline(infile,readline);
    
    //read in body,store all the information in bond_info as string type data
    while(getline(infile,readline)){
        stringstream ss(readline),new_ss;
        int Strnum=11;
        string str_data[Strnum];
        double coupon,price;
        
        //parser data //consider wrong data ,how will you handle it
        parser(ss,str_data,Strnum);
        new_ss<<str_data[2]<<'\t'<<str_data[9];
        new_ss>>coupon>>price;
        
        Maturity.push_back(str_data[3]);
        Coupon.push_back(coupon);
        Ask_price.push_back(price);
        Identifier.push_back(str_data[7]);
    }
    infile.close();
    
    int par_bond_num=int(Ask_price.size());
    vct_d Accrual(par_bond_num);
    vct_d time_to_maturity(par_bond_num);
    
    cout<<"Valuation_date "<<Valuation_date<<endl;
    
    time_to_Maturity(Maturity,Accrual,time_to_maturity,Valuation_date);
    
    for(int i=0;i<par_bond_num;i++){
        vector<double> c_t_p;//coupon,time to maturity,ask_price
        c_t_p.push_back(Coupon[i]);
        c_t_p.push_back(time_to_maturity[i]);
        c_t_p.push_back(Ask_price[i]);
        Bond_Info[Identifier[i]]=c_t_p;
    }

}

