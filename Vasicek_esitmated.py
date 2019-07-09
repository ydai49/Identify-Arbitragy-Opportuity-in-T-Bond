#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 21:26:59 2019

@author: daiyu
"""
import numpy as np
import scipy.optimize 


def Vasicek(theta,yeild):
    n=len(yeild)
    #print(theta)
    #print(yeild)
    b = (0,100)
    bnds = [b, b, b, b]
    def equation(x):
        #alpha,beta,sigma,spot
        res=np.zeros(n)
        infiniteR=x[1]-x[2]*x[2]/(2*x[0]*x[0])
        rest=0
        for i in range(n):
            res[i]=infiniteR-((infiniteR-x[3])*(1-np.exp(-x[0]*theta[i])))/(x[0]*theta[i]) \
            -x[2]*x[2]/(4*x[0]*x[0])*(1-np.exp(-theta[i]*x[0])*x[0])-yeild[i]
            rest=rest+res[i]*res[i]
        return(rest)
    result = scipy.optimize.minimize(equation, [yeild[1],yeild[n-1],yeild[2],yeild[0]],bounds=bnds,method='L-BFGS-B')
    
    print('Alpha: ',result.x[0])
    print('beta: ',result.x[1])
    print('sigma: ',result.x[2])
    print('short rate rt: ',result.x[3])
    
    return (list(result.x))

def Vasicek_mul(y2,y3,y4,y5,y6,y7,y8,y9,y10,theta):
    n=len(y2)#num of estiamteds 
    all_result=[]
    print(n)
    b = (0,100)
    bnds = [b, b, b, b]
    for i in range(n):
        yeild=[y2[i],y3[i],y4[i],y5[i],y6[i],y7[i],y8[i],y9[i],y10[i]]
        def equation(x):
            rest=0
            m=len(yeild)
            res=np.zeros(m)
            infiniteR=x[1]-x[2]*x[2]/(2*x[0]*x[0])
            for j in range(m):
                res[j]=infiniteR-((infiniteR-x[3])*(1-np.exp(-x[0]*theta[j])))/(x[0]*theta[j]) \
                -x[2]*x[2]/(4*x[0]*x[0])*(1-np.exp(-theta[j]*x[0])*x[0])-yeild[j]
                rest=rest+res[j]*res[j]
            return(rest)
        result = scipy.optimize.minimize(equation,[yeild[1],yeild[len(yeild)-1],yeild[2],yeild[0]],bounds=bnds,method='L-BFGS-B')
        all_result.append(result.x[0])
        all_result.append(result.x[1])
        all_result.append(result.x[2])
        all_result.append(result.x[3])
    np.array([0])
    np.array([0])
    print(result.x)
    
    return(list(all_result))
        
    
                