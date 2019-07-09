//
//  Vasicek_Estimation1.cpp
//  Vasicek_Model_Fitted_Yield
//
//  Created by 代雨 on 2019/5/2.
//  Copyright © 2019年 代雨. All rights reserved.
//

#include "Vasicek_Estimation1.hpp"
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

void datareader(string filename,map<double,double> &yield_data){
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
        yield_data[t]=y;
    }
    infile.close();
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

void Vasicek_estimation(map<double,double> yield_data,vector<double>& result,int num_of_result)
{
    // Initialize Python
    Py_Initialize();
    string path = "/Users/daiyu/Documents/2019 Spring Class/computational finance code";
    string chdir_cmd = string("sys.path.append(\"") + path + "\")";
    const char* cstr_cmd = chdir_cmd.c_str();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString(cstr_cmd);
    PyObject* moduleName = PyString_FromString("Vasicek_esitmated");
    PyObject* pModule = PyImport_Import(moduleName);
    // check if loading the module is succuessful
    if (!pModule)
    {
        cout << "[ERROR] Python get module failed." << endl;
        exit(1);
    }
    //cout << "[INFO] Python get module succeed." << endl;
    
    
    // Loading the function
    PyObject* pv = PyObject_GetAttrString(pModule, "Vasicek");
    if (!pv || !PyCallable_Check(pv))
    {
        cout << "[ERROR] Can't find funftion" << endl;
        exit(2);
    }
    cout << "[INFO] Get function (test_add) succeed." << endl;
    
    // pass the parameters
    // initialize a list
    PyObject*pyParams = PyList_New(0);
    PyObject*pyParams1 = PyList_New(0);
    
    for(auto &e: yield_data){
        //add member to the list
        PyList_Append(pyParams,Py_BuildValue("d",e.first));//time to maturity
        PyList_Append(pyParams1,Py_BuildValue("d",e.second));//yield
    }
    PyObject* args = PyTuple_New(2);
    PyTuple_SetItem(args,0,pyParams);
    PyTuple_SetItem(args,1,pyParams1);
    
    // make function call from python
    PyObject* pRet = PyObject_CallObject(pv, args);
    //cout<<"size:"<<size<<endl;
    // analysis return value from python
    if(pRet){
        for(int i = 0; i< num_of_result; ++i)
        {
            PyObject *pRet_new = PyList_GetItem(pRet, i);
            PyArg_Parse(pRet_new, "d", &result[i]);
        }
    }
    Py_Finalize();
}

void Vasicek_fitted_Yield_Curve(vector<double> result,map<double,double>& yield)
{
    double spot_rate=0.;
    double a= result[0];
    double b= result[1];
    double c= result[2];
    double r= result[3];
    int n=365*30;
    double dt=1.0/365;
    for(int i=1;i<=n;i++){
        double infinite=b-c*c/(2*a*a);
        double part1=(infinite-r)*(1-exp(-a*dt*i));
        double part2=c*c/(4*a*a);
        double part3=(1-exp(-a*dt*i))*(1-exp(-a*dt*i));
        spot_rate=infinite-1/(a*dt*i)*(part1-part2*part3);
        yield[i*dt]=spot_rate;
    }
}

void output_fitted_Yield_Curve(string filename,map<double,double> yield){
    ofstream outfile;
    outfile.open(filename);
    if(outfile.fail()){
        cerr<<"Open file fail"<<endl;
        exit(1);
    }
    outfile<<"time to maturity"<<'\t'<<"yield"<<endl;
    for(auto &e:yield)
        outfile<<e.first<<'\t'<<e.second<<endl;
    outfile.close();
}

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
void Vasicek_Bond_Price_1(map<string,vector<double>> Bond_Info,
                          map<string,vector<double>> &NewPrice,
                          vector<double> Vasicek_parameters,string filename)
{
    //Bond_Info: Ticker: coupon,time to maturity, price
    double a= Vasicek_parameters[0];
    double b= Vasicek_parameters[1];
    double c= Vasicek_parameters[2];
    double r= Vasicek_parameters[3];
    int num_bond=(int)Bond_Info.size();
    vector<string> keys;
    for(auto& e:Bond_Info)
        keys.push_back(e.first);
    for(int i=0;i<num_bond;i++){
        vector<double> new_old_price(2);
        //coupon pmt times
        new_old_price[0]=0;
        double maturity=Bond_Info[keys[i]][1];
        double hy=0.5;
        int coupon_times= floor(maturity/hy);
        double next_coupon_day=maturity-coupon_times*hy;
        double ai=(hy-next_coupon_day)/hy;//aacural
        vector<double> time_to_pmt;
        for(int j=0;j<=coupon_times;j++){
            time_to_pmt.push_back(next_coupon_day+j*hy);
        }
        vector<double>spot_rate(coupon_times+1);
        for(int j=0;j<=coupon_times;j++){
            double infinite=b-c*c/(2*a*a);
            double part1=(infinite-r)*(1-exp(-a*time_to_pmt[j]));
            double part2=c*c/(4*a*a);
            double part3=(1-exp(-a*time_to_pmt[j]))*(1-exp(-a*time_to_pmt[j]));
            spot_rate[j]=infinite-1/(a*time_to_pmt[j])*(part1-part2*part3);
            new_old_price[0]+=(Bond_Info[keys[i]][0]/2)/pow((1+spot_rate[j]),time_to_pmt[j]);
        }
        new_old_price[0]+=100./pow((1+spot_rate[coupon_times]),time_to_pmt[coupon_times]);
        new_old_price[0]=new_old_price[0]-ai*Bond_Info[keys[i]][0]/2;
        new_old_price[1]=Bond_Info[keys[i]][2];
        NewPrice[keys[i]]=new_old_price;
    }
    ofstream outfile;
    outfile.open(filename);
    outfile<<"Ticker"<<'\t'<<"Vasicek_price"<<'\t'<<"Mkt_Price"<<'\t'<<"DIff"<<endl;
    for(auto &a:NewPrice)
        if(abs(a.second[0]-a.second[1])>2){
            outfile<<a.first<<'\t'<<a.second[0]<<'\t'
            <<a.second[1]<<'\t'<<a.second[0]-a.second[1]<<endl;
        }
    outfile.close();
}
//**************************************************************
//this is a function that use Box-Muller Method to generate random number
double BoxMullerMethod(){
    double result;
    double x;
    double y;
    
    double xysquare;
    do{
        x=2.0*rand()/static_cast<double>(RAND_MAX)-1;
        y=2.0*rand()/static_cast<double>(RAND_MAX)-1;
        xysquare=x*x+y*y;
    }while(xysquare>=1.0);
    result =x*sqrt(-2*log(xysquare)/xysquare);
    return result;
}

void Vasicek_Simu1(vector<double> Vasicek_parameters,int num,int path,double dt,string filename)
{
    //Vasicek_parameters parameters: alpha,beta,sigma,rt
    double error;
    double alpha=Vasicek_parameters[0],beta=Vasicek_parameters[1];
    double sigma=Vasicek_parameters[2];
    double rate[path][num+1];
    double yield[path][num+1];

    for(int i=0;i<path;i++){
        rate[i][0]=Vasicek_parameters[3];
        for(int j=1;j<num+1;j++){
            error=BoxMullerMethod();
            rate[i][j]=rate[i][j-1]+(1-exp(-alpha*(dt*j)))/(alpha*(dt*j))*
                (alpha*(beta-rate[i][j-1])*dt+sigma*sqrt(dt)*error);
        }
    }
    ofstream outfile;
    outfile.open(filename);
    if(outfile.fail()){
        cerr<<"cannot open file"<<endl;
        exit(2);
    }
    for(int i=0;i<path;i++)
        outfile<<"Time"<<"Interest "<<i<<'\t';
    outfile<<endl;
    for(int j=0;j<num+1;j++){
        outfile<<dt*(j+1)<<'\t';
        for(int i=0;i<path;i++){
            outfile<<rate[i][j]<<'\t';
        }
        outfile<<endl;
    }
}



