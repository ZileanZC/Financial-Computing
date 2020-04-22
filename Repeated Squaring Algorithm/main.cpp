//
//  main.cpp
//  Repeated Squaring Algorithm
//
//  Created by Zilean Chen on 10/31/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
//
//Reference from Prof. R.S. Sreenivas's "strassen.h" for IE523: Financial Computing

#include <iostream>
#include <cmath>
#include <time.h>
#include <vector>
#include "newmat.h"
#include "newmatap.h"
#include "repeat_squaring.h"
using namespace std;

int main(int argc, const char * argv[]) {

    matrix_multiplication x;
    int exponent, num;
    sscanf (argv[1], "%d", &exponent);
    sscanf (argv[2], "%d", &num);
    //output what has been read
    cout<<"The number of rows/columns in the squared matrix is: "<<num<<endl;
    cout<<"The exponent is: "<<exponent<<endl;
    
    clock_t start_1,start_2,end_1,end_2;
    Matrix A=x.create_random_matrix(num);
    Matrix BFM_matrix, RSA_matrix;
    
    cout<<endl<<"When using Direct Multiplication:"<<endl;
    start_1=clock();
    BFM_matrix=x.Brute_force_multiply(A, exponent, num);
    x.print_matrix(BFM_matrix, num);
    end_1=clock();
    cout<<"This program took me "<<(double)(end_1-start_1)/CLOCKS_PER_SEC<<" seconds"<<endl;
    
    cout<<endl<<"When using Repeat Squaring Algorithm:"<<endl;
    start_2=clock();
    RSA_matrix=x.Repeat_squaring_algorithm(A, exponent, num);
    x.print_matrix(RSA_matrix, num);
    end_2=clock();
    cout<<"This program took me "<<(double)(end_2-start_2)/CLOCKS_PER_SEC<<" seconds"<<endl;
    
    //Double check result
    cout << "-------------------------------------------------------"<<endl;
    cout<<"To certify if they are identical, the sum differece is:"<<x.Judgement(BFM_matrix,RSA_matrix,num)<<endl;
    cout<<"The relative error of the multiplication is:"<<x.Judgement(BFM_matrix,RSA_matrix,num)/BFM_matrix(1,1)<<endl;
    
    //*
    vector <double> squaring_time;
    vector <double> direct_time;
    for (int i=0;i<300;i++)
    {
        clock_t start_1,start_2,end_1,end_2;
        start_1=clock();
        x.Brute_force_multiply(A, i, num);
        end_1=clock();
        direct_time.push_back((double)(end_1-start_1)/CLOCKS_PER_SEC);
        start_2=clock();
        x.Repeat_squaring_algorithm(A, i, num);
        end_2=clock();
        squaring_time.push_back((double)(end_2-start_2)/CLOCKS_PER_SEC);
    }
    ofstream s_t;
    ofstream d_t;
    s_t.open("squaring_time");
    d_t.open("direct_time");
    for (int i = 0; i <300; i++) {
        //cout<<direct_time[i]<<" ";
        s_t << squaring_time[i] << endl;
        d_t << direct_time[i] << endl;
    }
    //*/
}
