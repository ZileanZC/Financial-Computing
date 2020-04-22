//
//  main.cpp
//  Alice_Bob
//
//  Created by Zilean Chen on 10/25/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
//
// IE523: Financial Computation
// "How to lose as little as possible" by Addona, Wagon and Wilf
// Reference from 'hint.cpp' written by Prof. Sreenivas
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include "Alice_and_bob.h"
using namespace std;

#include <iostream>

int main(int argc, const char * argv[]) {
    // IE523: Financial Computation
    // "How to lose as little as possible" by Addona, Wagon and Wilf
    // Written by Prof. Sreenivas
    //
    I_have_nothing_apropos_for_this_class x;
    double alice_success_prob, bob_success_prob;
    sscanf (argv[1], "%lf", &alice_success_prob);
    sscanf (argv[2], "%lf", &bob_success_prob);
    
    cout << "Probability of success for Alice = " << alice_success_prob << endl;
    cout << "Probability of success for Bob = " << bob_success_prob << endl;
    
    x.set_probability(alice_success_prob, bob_success_prob);
    // I took time.h function reference from
    //https://blog.csdn.net/qq_36667170/article/details/79507547
    
    clock_t start,end;
    start=clock();
    
    int optimal = x.search_result();
    end=clock();
    /*long a,b;
    a=1e2;
    b=10e2;
    cout<<a<<"  "<<b<<endl;*/
    if (optimal > 0)
        cout << "The optimal number of coin tosses in each game is " << optimal << endl;
    else {
        cout << "The optimal number of coin tosses in each game exceeds 100... Quitting" << endl;
    }
    cout<<"This program took me "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;
    //data record
    x.create_data_file();
}

