//
//  main.cpp
//  Discrete Barrier
//
//  Created by Zilean Chen on 12/15/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
//

// Retrieved from: Prof. Ramavarapu
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
int no_of_trials, no_of_barriers;
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

double probability_pd (double final_stock_price){
    double a = initial_stock_price; double b = final_stock_price; double pd=1.0;
    for (int i=1;i<=no_of_barriers;i++){
        double mu=a+(i-1)*(b-a)/(no_of_barriers-1);
        double v=(i-1)*(no_of_barriers-i)/(no_of_barriers-1);
        double sd=sqrt(v);
        pd=pd*(1-N((barrier_price-mu)/sd));
    }
    return pd;
}

void explicit_simulation(){
    for (int i = 0; i < no_of_trials; i++)
    {
        double delta_T = expiration_time/((double)no_of_barriers);
        double delta_R = (risk_free_rate-0.5*pow(volatility,2))*delta_T;
        double delta_SD = volatility * sqrt(delta_T);
        
        int price_1_breach,price_2_breach,price_3_breach,price_4_breach;
        int t_1, t_2, t_3, t_4;
        double pd_1, pd_2, pd_3, pd_4;
        price_1_breach = 0; price_2_breach = 0; price_3_breach = 0; price_4_breach = 0;
        t_1=t_2=t_3=t_4=1;
        
        float S1 = initial_stock_price;
        float S2 = initial_stock_price;
        float S3 = initial_stock_price;
        float S4 = initial_stock_price;
        
        for (int j = 0; j < no_of_barriers; j++) {
            
            // create the unit normal variates using the Box-Muller Transform
            float x = get_uniform();
            float y = get_uniform();
            float a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            float b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            if (S1 <= barrier_price)
            {
                price_1_breach = 1; t_1 =0;
            }
            if (S2 <= barrier_price)
            {
                price_2_breach = 1; t_2 =0;
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
        pd_1 = probability_pd(S1); pd_2 = probability_pd(S2); pd_3 = probability_pd(S3); pd_4 = probability_pd(S4);
        explicit_call_price += (max(0.0, S1 - strike_price) * t_1 +max(0.0, S2 - strike_price) * t_2 + max(0.0, S3 - strike_price) * t_3 + max(0.0, S4 - strike_price) * t_4)/4.0;
        explicit_put_price += (max(0.0, -S1 + strike_price) * t_1 +max(0.0, -S2 + strike_price) * t_2 + max(0.0, -S3 + strike_price) * t_3 + max(0.0, -S4 + strike_price) * t_4)/4.0;
        
        
        adjust_call_price += (max(0.0, S1 - strike_price) * t_1 * pd_1 +max(0.0, S2 - strike_price) * t_2 * pd_2 + max(0.0, S3 - strike_price) * t_3 * pd_3 + max(0.0, S4 - strike_price) * t_4 * pd_4 )/4.0;
        adjust_put_price += (max(0.0, -S1 + strike_price) * t_1 * pd_1 +max(0.0, -S2 + strike_price) * t_2 * pd_2 + max(0.0, -S3 + strike_price) * t_3 * pd_3 + max(0.0, -S4 + strike_price) * t_4 * pd_4)/4.0;
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
    sscanf (argv[7], "%d", &no_of_barriers);
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
    cout << "Number of Discrete Barriers = "<<no_of_barriers <<endl;
    cout << "Barrier Price = "<<barrier_price<<endl;
    cout << "--------------------------------" << endl;
    cout << "--------------------------------" << endl;

    // Simulation
    explicit_simulation();
    // Print
    cout << "The average Call Price via explicit simulation of price paths = " << explicit_call_price << endl;
    cout << "The call price with Brownian-Bridge correction on the final price = " << adjust_call_price << endl;
    cout << "--------------------------------" << endl;
    cout << "The average Put Price via explicit simulation of price paths = " << explicit_put_price << endl;
    cout << "The put price with Brownian-Bridge correction on the final price = " << adjust_put_price << endl;
    cout << "--------------------------------" << endl;
}
