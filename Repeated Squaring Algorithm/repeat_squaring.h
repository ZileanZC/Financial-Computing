//
//  repeat_squaring.h
//  Repeated Squaring Algorithm
//
//  Created by Zilean Chen on 11/1/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
//

#ifndef repeat_squaring_h
#define repeat_squaring_h

#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include "newmat.h"
#include "newmatap.h"

using namespace std;

class matrix_multiplication
{
public:
    double get_uniform()
    {
        double tmp=(((double) random())/(pow(2.0, 31.0)-1.0));
        //[0,1]->(0,1)
        if (tmp==0||tmp==1)
        {tmp=get_uniform();}
        return tmp;
    }
    Matrix create_random_matrix(int n)
    {
        Matrix A(n,n);
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                double s=get_uniform();
                // 10*(0,1)-5->(5,5)
                A(i+1,j+1)=10.0*s-5.0;
            }
        }
        print_matrix(A, n);
        return A;
    }
    //print function modified from Prof.R.S. Sreenivas's "strassen.h"
    void print_matrix(Matrix A,int n)
    {
        cout << "The matrix is as follows:"<<std::endl;
        cout << "-------------------------------------------------------";
        cout<<endl;
        
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                cout<<setprecision(4)<<setw(12)<<A(i+1,j+1);
                //cout<<A(i+1,j+1);
                
            }
            cout<<endl;
        }
         
        cout << "-------------------------------------------------------";
        cout<<endl;
    }
    Matrix Brute_force_multiply(Matrix A, int e, int n)
    {
        Matrix BFM_matrix(n,n);
        BFM_matrix=A;
        for(int i=1;i<e;i++)
        {
            BFM_matrix=BFM_matrix*A;
        }
        return BFM_matrix;
    }
    Matrix Repeat_squaring_algorithm(Matrix A, int e, int n)
    {
        Matrix RSA_matrix(n,n);
        // Create diagonal matrix
        for (int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                if(i==j)
                    RSA_matrix(i+1,j+1)=1;
                else
                    RSA_matrix(i+1,j+1)=0;
            }
        }
        //print_matrix(RSA_matrix, n);
        // When e=0, return diagonal matrix
        if(e==0)
        {
            return RSA_matrix;
        }
        if (e%2==1)
        {
            return (A*Repeat_squaring_algorithm(A*A, (e-1)/2, n));
        }
        else
        {
            return (Repeat_squaring_algorithm(A*A, e/2, n));
        }
    }
    double Judgement(Matrix M, Matrix N, int n)
    {
        double err=0;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                err+=M(i+1,j+1)-N(i+1,j+1);
            }
        }
        return err;
    }

};

#endif /* repeat_squaring_h */
