//
//  Alice_and_bob.h
//  Alice_Bob
//
//  Created by Zilean Chen on 10/25/19.
//  Copyright © 2019 Zilean Chen. All rights reserved.
// Reference from 'hint.h' written by Prof. Sreenivas
//

#ifndef Alice_and_bob_h
#define Alice_and_bob_h

#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

class I_have_nothing_apropos_for_this_class
{
private:
    double alice_probability, bob_probability;
    
    // private member function: uniform RV generator
    double get_uniform()
    {
        return (((double) random())/(pow(2.0, 31.0)-1.0));
    }
    
    // private member function: nCi (i.e. n-take-i)
    long int take(int n, int i)
    {
        // Get rid of conception error
        if(i>n)
        {
            cout<<"Wrong input in taking "<<n<<" from "<<i<<"."<<endl;
            return 0;
        }
        // write a **RECURSIVE** implementation of n-take-i.
        // If you made it non-recurisive (i.e. n!/((n-i)!i!)) -- it
        // will take too long for large sizes
        else
        {
            // take(n,n)=1 && take(n,0)=1
            if(n==i || i==0)
            {
                return 1;
            }
            // take(n,i)=take(n-1,i)+take(n-1,i-1)
            else
            {
                return take(n-1,i)+take(n-1,i-1);
            }
        }
    }
    
    // this routine implements the probability that Alice has more
    // heads than Bob after n-many coin tosses
    double theoretical_value(double q, double p, int n)
    {
        // implement equation 1.1 of Addona-Wagon-Wilf paper
        double Bob_head=0;
        for (int r=0;r<n;r++)
        {
            double Alice_head=0;
            for (int s=r+1;s<=n;s++)
            {
                Alice_head=Alice_head+take(n,s)*pow(q,s)*pow((1-q),n-s);
            }
            Bob_head=Bob_head+take(n,r)*pow(p,r)*pow((1-p),n-r)*Alice_head;
        }
        double Alice_win=Bob_head;
        return Alice_win;
    }

public:
    // public function:
    void set_probability(double alice_p, double bob_p)
    {
        alice_probability = alice_p;
        bob_probability = bob_p;
    }
    
    // probability of Alice winning the game.
    double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
    {
        int no_of_wins_for_alice = 0;
        for (int i = 0; i < no_of_trials; i++)
        {
            int number_of_heads_for_alice = 0;
            int number_of_heads_for_bob = 0;
            for (int j = 0; j < number_of_coin_tosses_in_each_game; j++)
            {
                if (get_uniform() < alice_probability)
                    number_of_heads_for_alice++;
                if (get_uniform() < bob_probability)
                    number_of_heads_for_bob++;
            }
            if (number_of_heads_for_alice > number_of_heads_for_bob)
                no_of_wins_for_alice++;
        }
        return (((double) no_of_wins_for_alice)/((double) no_of_trials));
    }
        
    //data record
    vector <double> simulated_data;
    vector <double> theoretical_data;
    
    void create_data_file()
    {
        for (int i = 0; i < 30; i++)
        {
            theoretical_data.push_back(theoretical_value(alice_probability, bob_probability, i));
            simulated_data.push_back(simulated_value(i, 1e6)); //here modify the no_of_trials
        }
           
        ofstream simulated_file;
        ofstream theoretical_file;
        simulated_file.open("simulated_data");
        theoretical_file.open("theoretical_data");
        
        long int N=simulated_data.size();
        for (int i = 0; i < N; i++)
        {
            simulated_file << simulated_data[i] << endl;
            theoretical_file << theoretical_data[i] << endl;
        }
       }
    
    //optimization result
    int search_result()
    {
        // implememt a discrete-search procedure for the optimal n-value.
        // start with n = 1 and find the discrete-value of n that has
        // the largest probability for Alice winning.  Why would this work?
        // See Theorem 2.2 of the paper for the reason!
        for (int i=0;i<=1e2;i++)
        {
            // f(n)>=f(n-1) && f(n)>=f(n+1)
            if ((theoretical_value(alice_probability, bob_probability, i) >= theoretical_value(alice_probability, bob_probability, i - 1)) &&
            (theoretical_value(alice_probability, bob_probability, i) >= theoretical_value(alice_probability, bob_probability, i + 1)))
                // In this case, we need to find out the smallest n
                // If we want to find out all n(s), we can assign an int *array=new int [n]；
                return i;
        }
        return 0;
    }
};


#endif /* Alice_and_bob_h */






