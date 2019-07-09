//
//  Nelson Seigel.cpp
//  Grid Search For Nelson Seigel and yield Fit
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#include "Nelson Seigel.hpp"
#include <fstream>
#include <iostream>
#include<string>
#include<math.h>
#include<vector>
#include<map>
#include<boost/date_time/gregorian/gregorian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
using namespace std;
using vct_d=vector<double>;
using d_matrix=vector<vct_d>;

void datareader(string filename,map<double,double> &yield_data,map<double,double> &all_data){
    ifstream infile;
    infile.open(filename);
    if(infile.fail()){
        cerr<<"Can't open file"<<endl;
        exit(2);
    }
    string readline;
    getline(infile,readline);
    cout<<"Data chosen to do estimation"<<endl;
    cout<<readline<<endl;
    
    vector<double> time;
    vector<double> yield;
    
    while(getline(infile,readline)){
        stringstream ss(readline),ssr;
        double t,y;
        parser(ss, t, y);
        time.push_back(t);
        yield.push_back(y);
        all_data[t]=y;
    }
    infile.close();
    // find target data:
    int n=(int)time.size();
    for(int i=0;i<n;i++)
        if(time[i]>0.5){
            yield_data[time[i]]=yield[i];
            break;}
    for(int i=0;i<n;i++)
        if(time[i]>1.0){
            yield_data[time[i]]=yield[i];
            break;}
    for(int i=0;i<n;i++)
        if(time[i]>2.0){
            yield_data[time[i]]=yield[i];
            break;}
    for(int i=0;i<n;i++)
        if(time[i]>3.0){
            yield_data[time[i]]=yield[i];
            break;}
    for(int i=0;i<n;i++)
        if(time[i]>5.0){
            yield_data[time[i]]=yield[i];
            break;}
    for(int i=0;i<n;i++)
        if(time[i]>7.0){
            yield_data[time[i]]=yield[i];
            break;}
    for(int i=0;i<n;i++)
        if(time[i]>10.0){
            yield_data[time[i]]=yield[i];
            break;}
    for(auto & e:yield_data)
        cout<<e.first<<' ' <<e.second<<endl;
    cout<<"Total number of choosing data: "<<yield_data.size()<<endl;
}

void parser(stringstream& ss,double &time,double &yield){
    string t,y;
    getline(ss,t,',');
    getline(ss,y,',');
    stringstream ssr;
    ssr<<t<<'\t'<<y;
    ssr>>time>>yield;
}

//step 2: Grid Search to find the proper lambda for NS
void Grid_Search_NS(map<double,double> yield_data,map<double,double> all_data,vector<double> &parameters){
//
    double lambda=0.6;
    vector<double> grid;
    double low=0;
    double up=30;
    double dg=0.001;
    int n=(up-low)/dg;
    for(int i=0;i<n;i++)
        grid.push_back(i*dg);
    int n1=(int)yield_data.size();
    
    vector<double>time;
    for(auto&e:yield_data)
        time.push_back(e.first);
    
    boost::numeric::ublas::matrix<double> X(n1,3);
    boost::numeric::ublas::vector<double> Y(n1);
    boost::numeric::ublas::vector<double> new_Y(n1);
    boost::numeric::ublas::vector<double> beta(3);
    boost::numeric::ublas::matrix<double> all_beta_grid(n,3);
    //filling X with factor_loading data
    for(int k=0;k<n;k++){
        lambda=grid[k];
        for(int i=0;i<n1;i++){
            X(i,0)=1.0;
            X(i,1)=(1-exp(-lambda*time[i]))/(lambda*time[i]);
            X(i,2)=(1-exp(-lambda*time[i]))/(lambda*time[i])-
            exp(-lambda*time[i]);
            Y(i)=yield_data[time[i]];
        }
        boost::numeric::ublas::matrix<double> trans_X;
        boost::numeric::ublas::matrix<double> mul_XX;
        trans_X=trans(X);
        mul_XX=prod(trans_X,X);
        new_Y=prod(trans_X,Y);
        boost::numeric::ublas::permutation_matrix<double> P(3);
        lu_factorize(mul_XX,P);
        // Now mul_XX and P contain the LU factorization of mul_XX
        beta = new_Y;
        lu_substitute(mul_XX,P,beta);
        all_beta_grid(k,0)=beta(0);
        all_beta_grid(k,1)=beta(1);
        all_beta_grid(k,2)=beta(2);
    }
    // calculated mse for each parameters
    int n_all=(int)all_data.size();
    boost::numeric::ublas::matrix<double> mse_X(n_all,3);
    boost::numeric::ublas::vector<double> mse_Y(n_all);
    boost::numeric::ublas::vector<double> mse_beta(3);
    vector<double> mse(n);
    vector<double> mse_time;
    for(auto &e:all_data)
        mse_time.push_back(e.first);
    
    for(int i=0;i<n;i++){
        lambda=grid[i];
        //filling the factor loading matrix
        for(int k=0;k<n_all;k++){
            mse_X(k,0)=1.0;
            mse_X(k,1)=(1-exp(-lambda*mse_time[k]))/(lambda*mse_time[k]);
            mse_X(k,2)=(1-exp(-lambda*mse_time[k]))/(lambda*mse_time[k])-
            exp(-lambda*mse_time[k]);
            mse_Y(k)=all_data[mse_time[k]];//true yield
        }
        //cout<<mse_Y<<endl;
        double sum=0;
        //filling the factor with estimated value
        mse_beta(0)=all_beta_grid(i,0);
        mse_beta(1)=all_beta_grid(i,1);
        mse_beta(2)=all_beta_grid(i,2);
        boost::numeric::ublas::vector<double> estimated_Y;
        estimated_Y=prod(mse_X,mse_beta);
        for(int j=0;j<n_all;j++)
            sum+=(estimated_Y(j)-mse_Y(j))*(estimated_Y(j)-mse_Y(j));
        mse[i]=sum/n_all;
    }
    double min=mse[1];
    int minIndex=0;
    for(int i=1;i<(int)mse.size();i++){
        if(min>mse[i]){
            min=mse[i];
            minIndex=i;
        }
    }
    cout<<"In the grid search to find best lambda:\nThe upper bound: "<<up<<"\nThe lower bound: 0\n";
    lambda=grid[minIndex];
    cout<<"The best fit lambda is: "<<lambda<<"\nThe MSE is: "<<mse[minIndex]<<endl;
    
    //reestimated beta with best fit lambda and all the data
    for(int k=0;k<n_all;k++){
        mse_X(k,0)=1.0;
        mse_X(k,1)=(1-exp(-lambda*mse_time[k]))/(lambda*mse_time[k]);
        mse_X(k,2)=(1-exp(-lambda*mse_time[k]))/(lambda*mse_time[k])-
        exp(-lambda*mse_time[k]);
        mse_Y(k)=all_data[mse_time[k]];//true yield
    }
    boost::numeric::ublas::permutation_matrix<double> P(3);
    boost::numeric::ublas::matrix<double> trans_all_X;
    boost::numeric::ublas::matrix<double> mul_all_XX;
    boost::numeric::ublas::vector<double> all_Y;
    trans_all_X=trans(mse_X);
    
    mul_all_XX=prod(trans_all_X,mse_X);
    
    all_Y=prod(trans_all_X,mse_Y);
    
    lu_factorize(mul_all_XX,P);
    // Now mul_XX and P contain the LU factorization of mul_XX
    lu_substitute(mul_all_XX,P,all_Y);
    
    parameters[0]=all_Y(0);
    parameters[1]=all_Y(1);
    parameters[2]=all_Y(2);
    parameters[3]=lambda;
    cout<<"beta0: "<<parameters[0]
    <<"\nbeta1: "<<parameters[1]
    <<"\nbeta2: "<<parameters[2]<<endl;
}

void NS_fit_yield_curve(vector<double> &parameters,map<double,double> &yield){
    int n=365*30;
    double dt=1.0/365;
    boost::numeric::ublas::matrix<double> X(n,3);
    boost::numeric::ublas::vector<double> Y(n);
    boost::numeric::ublas::vector<double> beta(3);
    beta(0)=parameters[0];
    beta(1)=parameters[1];
    beta(2)=parameters[2];
    double lambda=parameters[3];
    for(int i=0;i<n;i++){
        X(i,0)=1.0;
        X(i,1)=(1-exp(-lambda*(1+i)*dt))/(lambda*(1+i)*dt);
        X(i,2)=(1-exp(-lambda*(1+i)*dt))/(lambda*(1+i)*dt)-
            exp(-lambda*(1+i)*dt);
        }
    Y=prod(X, beta);
    for(int i=0;i<n;i++)
        yield[(i+1)*dt]=Y(i);
}

void dataoutput(string filename,map<double,double> &yield){
    ofstream outfile;
    outfile.open(filename);
    if(outfile.fail()){
        cerr<<"can't not open file"<<endl;
        exit(1);
    }
    outfile<<"time to maturity"<<'\t'<<"yield"<<endl;
    for(auto & e:yield)
        outfile<<e.first<<'\t'<<e.second<<endl;
    outfile.close();
}

//Step 5
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
        parser_bond(ss,str_data,Strnum);
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

void parser_bond(stringstream & ss,string * str_data,int str_num)
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
    int Diff_Day[int(Maturity.size())];
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
//Step 5:

void N_S_Bond_Price(map<string,vector<double>> Bond_Info,
                    map<string,vector<double>>& NewPrice,vector<double> NS_parameters,string filename)
{
    //Bond_Info: Ticker: coupon,time to maturity, price
    int num_bond=(int)Bond_Info.size();
    vector<string> keys;
    for(auto& e:Bond_Info)
        keys.push_back(e.first);
    boost::numeric::ublas::vector<double> beta(3);
    beta(0)=NS_parameters[0];
    beta(1)=NS_parameters[1];
    beta(2)=NS_parameters[2];
    double lambda=NS_parameters[3];
    for(int i=0;i<num_bond;i++){
        double maturity=Bond_Info[keys[i]][1];
        double hy=0.5;
        int coupon_times= floor(maturity/hy);
        double next_coupon_day=maturity-coupon_times*hy;
        double ai=(hy-next_coupon_day)/hy;//accural
        vector<double> time_to_pmt;
        for(int j=0;j<=coupon_times;j++){
            time_to_pmt.push_back(next_coupon_day+j*hy);
        }
        boost::numeric::ublas::matrix<double> X(coupon_times+1,3);
        boost::numeric::ublas::vector<double> Y(coupon_times+1);
        boost::numeric::ublas::vector<double> spot(coupon_times+1);
        for(int j=0;j<coupon_times+1;j++){
            X(j,0)=1.0;
            X(j,1)=(1-exp(-lambda*time_to_pmt[j]))/(lambda*time_to_pmt[j]);
            //cout<<(1-exp(-lambda*time_to_pmt[j]))<<endl;
            X(j,2)=(1-exp(-lambda*time_to_pmt[j]))/(lambda*time_to_pmt[j])-
            exp(-lambda*time_to_pmt[j]);
        }
        
        Y=prod(X, beta);
        //Dicount each coupon back to PV
        double dis_price=0;
        for(int j=0;j<coupon_times+1;j++)
            dis_price+=(Bond_Info[keys[i]][0]/2)/pow((1+Y(j)),time_to_pmt[j]);
        dis_price+=100.0/pow((1+Y(coupon_times)),time_to_pmt[coupon_times]);
        dis_price=dis_price-ai*Bond_Info[keys[i]][0]/2;
        vector<double> np;
        np.push_back(dis_price);
        np.push_back(Bond_Info[keys[i]][2]);
        NewPrice[keys[i]]=np;
    }
    ofstream outfile;
    outfile.open(filename);
    outfile<<"Ticker"<<'\t'<<"NS_price"<<'\t'<<"Mkt_Price"<<'\t'<<"DIff"<<endl;
    for(auto &a:NewPrice)
        if(abs(a.second[0]-a.second[1])>2){
            outfile<<a.first<<'\t'<<a.second[0]<<'\t'
            <<a.second[1]<<'\t'<<a.second[0]-a.second[1]<<endl;
        }
    outfile.close();
}

