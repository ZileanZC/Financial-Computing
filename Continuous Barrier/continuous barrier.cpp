//
//  main.cpp
//  Continuous Barrier
//
//  Created by Zilean Chen on 12/15/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
//

// Retrieved from : Prof. Sreenivas

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
double barrier_price;
int no_of_trials, no_of_points;
double explicit_call_price, explicit_put_price, adjust_call_price, adjust_put_price;

// Cumulative Normal Distribution
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

double max(double a, double b) {
    return (b < a )? a:b;
}

// Denote to generate uniform distribution
double get_uniform()
{
    return (((double) random())/(pow(2.0, 31.0)-1.0));
}

// Black-Scholes Model Call/Put Option Price
double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

float closed_form_down_and_out_european_call_option()
{
    // I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
    float K = (2*risk_free_rate)/(volatility*volatility);
    float A = option_price_call_black_scholes(initial_stock_price, strike_price,
                                              risk_free_rate, volatility, expiration_time);
    float B = (barrier_price*barrier_price)/initial_stock_price;
    float C = pow(initial_stock_price/barrier_price, -(K-1));
    float D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
    return (A - D*C);
}

float closed_form_down_and_in_european_put_option()
{
    // just making it easier by renaming the global variables locally
    float S = initial_stock_price;
    float r = risk_free_rate;
    float T = expiration_time;
    float sigma = volatility;
    float H = barrier_price;
    float X = strike_price;
    
    // Took these formulae from some online reference
    float lambda = (r+((sigma*sigma)/2))/(sigma*sigma);
    float temp = 2*lambda - 2.0;
    float x1 = (log(S/H)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    float y = (log(H*H/(S*X))/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    float y1 = (log(H/S)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    return (-S*N(-x1) + X*exp(-r*T)*N(-x1 + sigma*sqrt(T)) +
            S*pow(H/S, 2*lambda)*(N(y)-N(y1)) -
            X*exp(-r*T)*pow(H/S, temp)*(N(y-sigma*sqrt(T))-N(y1-sigma*sqrt(T))));
}

float closed_form_down_and_out_european_put_option()
{
    float vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
                                                       risk_free_rate, volatility, expiration_time);
    float put_down_in = closed_form_down_and_in_european_put_option();
    return (vanilla_put - put_down_in);
}

// The probability that an asset breaches the barrier price

double probability_pc (double final_stock_price){

    double tmp = 0.0;
    if (initial_stock_price<= barrier_price || final_stock_price <= barrier_price)
        tmp=1.0;

    else
        tmp = exp(-2*log(initial_stock_price/barrier_price)*log(final_stock_price/barrier_price)/(pow(volatility,2)*expiration_time));
    return tmp;
}

void explicit_simulation(){
    for (int i = 0; i < no_of_trials; i++)
    {
        double delta_T = expiration_time/no_of_points;
        double delta_R = (risk_free_rate-0.5*pow(volatility,2))*delta_T;
        double delta_SD = volatility * sqrt(delta_T);
        
        int price_1_breach,price_2_breach,price_3_breach,price_4_breach;
        int t_1, t_2, t_3, t_4;
        double pc_1, pc_2, pc_3, pc_4;
        price_1_breach = 0; price_2_breach = 0; price_3_breach = 0; price_4_breach = 0;
        t_1=t_2=t_3=t_4=1;
        
        float S1 = initial_stock_price;
        float S2 = initial_stock_price;
        float S3 = initial_stock_price;
        float S4 = initial_stock_price;
        
        for (int j = 0; j < no_of_points; j++) {
            
            // create the unit normal variates using the Box-Muller Transform
            float x = get_uniform();
            float y = get_uniform();
            float a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            float b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            if (S1 <= barrier_price)
            {
                price_1_breach = 1; t_1 = 0;
            }
            if (S2 <= barrier_price)
            {
                price_2_breach = 1; t_2 = 0;
            }
            if (S3 <= barrier_price)
            {
                price_3_breach = 1; t_3=0;
            }
            if (S4 <= barrier_price)
            {
                price_4_breach = 1; t_4=0;
            }
            S1 = S1*exp(delta_R + delta_SD*a);
            S2 = S2*exp(delta_R - delta_SD*a);
            S3 = S3*exp(delta_R + delta_SD*b);
            S4 = S4*exp(delta_R - delta_SD*b);
        }
        pc_1 = probability_pc(S1); pc_2 = probability_pc(S2); pc_3 = probability_pc(S3); pc_4 = probability_pc(S4);
        explicit_call_price += (max(0.0, S1 - strike_price) * t_1 +max(0.0, S2 - strike_price) * t_2 + max(0.0, S3 - strike_price) * t_3 + max(0.0, S4 - strike_price) * t_4)/4.0;
        explicit_put_price += (max(0.0, -S1 + strike_price) * t_1 +max(0.0, -S2 + strike_price) * t_2 + max(0.0, -S3 + strike_price) * t_3 + max(0.0, -S4 + strike_price) * t_4)/4.0;
        
        
        adjust_call_price += (max(0.0, S1 - strike_price) * t_1 * (1-pc_1) +max(0.0, S2 - strike_price) * t_2 * (1-pc_2) + max(0.0, S3 - strike_price) * t_3 * (1-pc_3) + max(0.0, S4 - strike_price) * t_4 * (1-pc_4))/4.0;
        adjust_put_price += (max(0.0, -S1 + strike_price) * t_1 * (1-pc_1) +max(0.0, -S2 + strike_price) * t_2 * (1-pc_2) + max(0.0, -S3 + strike_price) * t_3 * (1-pc_3) + max(0.0, -S4 + strike_price) * t_4 *(1-pc_4))/4.0;
    }
    explicit_call_price = (exp(-risk_free_rate*expiration_time)*(explicit_call_price/((double)no_of_trials)));
    explicit_put_price = (exp(-risk_free_rate*expiration_time)*(explicit_put_price/((double)no_of_trials)));
    
    adjust_call_price = (exp(-risk_free_rate*expiration_time)*(adjust_call_price/((double)no_of_trials)));
    adjust_put_price = (exp(-risk_free_rate*expiration_time)*(adjust_put_price/((double)no_of_trials)));
}

int main (int argc, char* argv[])
{
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &no_of_trials);
    sscanf (argv[7], "%d", &no_of_points);
    sscanf (argv[8], "%lf", &barrier_price);
    
    double R = (risk_free_rate - 0.5*pow(volatility,2))*expiration_time;
    double SD = volatility*sqrt(expiration_time);
    
    cout << "--------------------------------" << endl;
    cout << "European Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Divisions = "<<no_of_points <<endl;
    cout << "Barrier Price = "<<barrier_price<<endl;
    cout << "--------------------------------" << endl;
    cout << "--------------------------------" << endl;

    // Simulation
    explicit_simulation();
    // Print
    cout << "The average Call Price by explicit simulation = " << explicit_call_price << endl;
    cout << "The call price using the (1-p)-adjustment term = " << adjust_call_price << endl;
    cout << "Theoretical Call Price = " << closed_form_down_and_out_european_call_option() << endl;
    cout << "--------------------------------" << endl;
    cout << "The average Put Price by explicit simulation = " << explicit_put_price << endl;
    cout << "The put price using the (1-p)_adjustment term = " << adjust_put_price << endl;
    cout << "Theoretical Put Price = " << closed_form_down_and_out_european_put_option() << endl;
    cout << "--------------------------------" << endl;
}
