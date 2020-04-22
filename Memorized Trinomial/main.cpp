//
//  main.cpp
//  Memorized Trinomial
//
//  Created by Zilean Chen on 11/14/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
//
//Reference comes from
//black_scholes.h
//  originally written by Prof. Sreenivas

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "newmat.h"
using namespace std;

double up_factor, uptick_prob, risk_free_rate, strike_price;
double downtick_prob, notick_prob;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

//Black-Scholes Model
double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
}
double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
}
double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
}

//Memozied European Option
double european_call_option(int k, int i, double** memoization) {
    if (k == no_of_divisions) {
        return max(0.0, (initial_stock_price * pow(up_factor, i) - strike_price));
    }
    if (memoization[k][no_of_divisions + i] == -1) {
        memoization[k][no_of_divisions + i] =
            (uptick_prob * european_call_option(k + 1, i + 1, memoization) +
                notick_prob * european_call_option(k + 1, i, memoization) +
                downtick_prob * european_call_option(k + 1, i - 1, memoization)) / R;
        return memoization[k][no_of_divisions + i];
    }
    else
        return memoization[k][no_of_divisions + i];
}

double european_put_option(int k, int i, double** memoization) {
    if (k == no_of_divisions) {
        return max(0.0, strike_price - (initial_stock_price * pow(up_factor, i)));
    }
    if (memoization[k][no_of_divisions + i] == -1) {
        memoization[k][no_of_divisions + i] =
            (uptick_prob * european_put_option(k + 1, i + 1, memoization) +
                notick_prob * european_put_option(k + 1, i, memoization) +
                downtick_prob * european_put_option(k + 1, i - 1, memoization)) / R;
        return memoization[k][no_of_divisions + i];
    }
    else
        return memoization[k][no_of_divisions + i];
}

//Memoized American Option

double american_call_option(int k, int i, double** memoization)
{
    if (k == no_of_divisions)
        return max(0.0, (initial_stock_price * pow(up_factor, i) - strike_price));
    else
    {
        if (memoization[k][no_of_divisions + i] == -1)
        {
            memoization[k][no_of_divisions + i] = max((initial_stock_price * pow(up_factor, i) - strike_price),
                (uptick_prob * american_call_option(k + 1, i + 1, memoization) +
                    downtick_prob * american_call_option(k + 1, i - 1, memoization) + notick_prob * american_call_option(k + 1, i, memoization)) / R);
            return memoization[k][no_of_divisions + i];
        }
        else
            return memoization[k][no_of_divisions + i];
    }
}

double american_put_option(int k, int i, double** memoization)
{
    if (k == no_of_divisions)
        return max(0.0, (strike_price - initial_stock_price * pow(up_factor, i)));
    else
    {
        if (memoization[k][no_of_divisions + i] == -1)
        {
            memoization[k][no_of_divisions + i] = max((strike_price - initial_stock_price * pow(up_factor, i)),
                (uptick_prob * american_put_option(k + 1, i + 1, memoization) +
                    downtick_prob * american_put_option(k + 1, i - 1, memoization) + notick_prob * american_put_option(k + 1, i, memoization)) / R);
            return memoization[k][no_of_divisions + i];
        }
        else
            return memoization[k][no_of_divisions + i];
    }
}



int main(int argc, char *argv[]){
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%lf", &risk_free_rate);
    sscanf (argv[4], "%lf", &volatility);
    sscanf (argv[5], "%lf", &initial_stock_price);
    sscanf (argv[6], "%lf", &strike_price);
    
    up_factor = exp(volatility*sqrt(2*expiration_time/((float) no_of_divisions)));
    R = exp(risk_free_rate*expiration_time/((float) no_of_divisions));
    uptick_prob = (sqrt(R) - (1/sqrt(up_factor)))/(sqrt(up_factor)-(1/sqrt(up_factor)));
    uptick_prob=pow(uptick_prob,2);
    downtick_prob=(sqrt(up_factor)-sqrt(R))/(sqrt(up_factor)-(1/sqrt(up_factor)));
    downtick_prob=pow(downtick_prob,2);
    notick_prob=1-uptick_prob-downtick_prob;
    
    //initialization
    
    int n=no_of_divisions;
    double ** memoization_1, ** memoization_2, ** memoization_3, ** memoization_4;
    memoization_1 = new double* [n];
    memoization_2 = new double* [n];
    memoization_3 = new double* [n];
    memoization_4 = new double* [n];
    for (int i = 0; i < n;i++)
    {
        memoization_1[i] = new double[2 * n + 1];
        memoization_2[i] = new double[2 * n + 1];
        memoization_3[i] = new double[2 * n + 1];
        memoization_4[i] = new double[2 * n + 1];
        
    }
    for (int i = 0; i < n ;i++)
    {
        for (int j = 0; j <2* n + 1; j++)
        {
            memoization_1[i][j] = -1;
            memoization_2[i][j] = -1;
            memoization_3[i][j] = -1;
            memoization_4[i][j] = -1;
        }
    }
    
    //European Option Price
    
    cout << "(Memoized) Recursive Trinomial European Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    //cout << "R = " << R << endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout<<"Downtick Probability = "<<downtick_prob<<endl;
    cout<<"Notick Probability = "<<notick_prob<<endl;
    cout << "--------------------------------------" << endl;
    
    //european call price
    double call_price = european_call_option(0, 0, memoization_1);
    cout << "Trinomial Price of an European Call Option = " << call_price << endl;
    double call_price_BS=option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,volatility, expiration_time);
    cout<<"Call Price according to Black-Scholes: "<<call_price_BS<<endl;
    cout << "--------------------------------------" << endl;
    //european put price
    double put_price = european_put_option(0, 0, memoization_2);
    cout << "Trinomial Price of an European Put Option = " << put_price << endl;
    double put_price_BS=option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,volatility, expiration_time);
    cout<<"Put Price according to Black-Scholes: "<<put_price_BS<<endl;
    cout << "--------------------------------------" << endl;
    //Put-Call Parity
    cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
    cout <<  initial_stock_price << " + " << put_price << " - " << call_price;
    cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
    cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
    if (abs(initial_stock_price + put_price - call_price - strike_price*exp(-risk_free_rate*expiration_time)) <= 1e-3)
    {cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
        cout << initial_stock_price + put_price - call_price << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;}
    else
        cout << "Looks like Put-Call Parity does NOT hold" << endl;
    cout << "--------------------------------------" << endl;
    
    //American Option Price

    cout << endl;
    cout << "(Memoized) Recursive Trinomial American Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    //cout << "R = " << R << endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    cout<<"Downtick Probability = "<<downtick_prob<<endl;
    cout<<"Notick Probability = "<<notick_prob<<endl;
    cout << "--------------------------------------" << endl;
    
    //american call option
    double call_price_US=american_call_option(0, 0, memoization_3);
    cout << "Trinomial Price of an American Call Option = " << call_price_US << endl;
    //american put option
    double put_price_US=american_put_option(0,0,memoization_4);
    cout << "Trinomial Price of an European Call Option = " << put_price_US << endl;
    cout << "--------------------------------------" << endl;
    
    //delete pointer-to-pointer

    delete [] memoization_1;
    delete [] memoization_2;
    delete [] memoization_3;
    delete [] memoization_4;

}
