//
//  main.cpp
//  RPMD with reduced units
//
//  Created by Harry Ryu on 2017. 5. 11..
//  Copyright © 2017년 Harry Ryu. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

int P;
// Number of beads
double h_bar = 1;
// In en*fs
double mass = 1;
// In au (1 au = 1.660539e-27 kg)
// Set it as the mass of proton
int N;
// Number of timesteps taken
// Make this even to make life easier
int corr_length;
double beta;
double beta_P = beta/P;
// In 1/en
double dt;
// Timestep
int N_samp;
// Number of sampling from the canonical distribution
double quad;
double cub;
double quart;

double momentum_dist()
// Box Muller algorithm that returns Gaussian random number with 0 mean and nonunitary variance
// Variance for the momentum in Boltzmann distirbution is mass/beta
// Return momentum in unit of amu*Å/fs
{
    double first, v1, v2, rsq, fac;
    int again = 1;
    while (again)
    {
        v1 = 2.0*rand()/RAND_MAX - 1.0;
        v2 = 2.0*rand()/RAND_MAX - 1.0;
        rsq = v1*v1 + v2*v2;
        if (rsq < 1.0 & rsq != 0.0) again = 0;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    first = v1 * fac * sqrt(mass/beta_P);
    return first;
}

double pot_deriv(double x)
{
    double value;
    value = 2 * quad * x + 3 * cub * x * x + 4 * quart * x * x * x;
    return value;
}

void transform(double *output, double **matrix, double *input)
{
    for (int i = 0 ; i < P ; ++i)
    {
        output[i] = 0.0;
        for (int j = 0 ; j < P ; ++j)
        {
            output[i] = output[i] + matrix[i][j] * input[j];
        }
    }
}

void bead_mom_update_half(double *mom_bead, double *pos_bead)
{
    for (int i = 0 ; i < P ; ++i)
    {
        mom_bead[i] = mom_bead[i] - (dt / 2) * pot_deriv(pos_bead[i]);
    }
}

void normal_update(double *mom_new, double *pos_new, double *mom_old, double *pos_old, double *freq)
{
    mom_new[0] = mom_old[0];
    pos_new[0] = (dt / mass) * mom_old[0] + pos_old[0];
    for (int i = 1 ; i < P ; ++i)
    {
        mom_new[i] = cos(freq[i] * dt) * mom_old[i] - mass * freq[i] * sin(freq[i] * dt) * pos_old[i];
        pos_new[i] = sin(freq[i] * dt) * mom_old[i] / (mass * freq[i]) + cos(freq[i] * dt) * pos_old[i];
    }
}

double *calc_corr(double *q)
{
    double *corr = new double[corr_length];
    for (int n = 0 ; n < corr_length ; ++n)
    {
        corr[n] = 0.0;
        for (int i = 0 ; i < (N-n) ; ++i)
        {
            corr[n] = corr[n] + q[i] * q[i+n] / (N - n);
        }
    }
    return corr;
}

int main()
{
    char str;
    
    FILE * pFile;
    
    pFile = fopen ("input.txt", "r"); // Reading in the parameters
    
    fscanf(pFile, "%i\n", &P);
    fscanf(pFile, "%i\n", &N);
    fscanf(pFile, "%i\n", &corr_length);
    fscanf(pFile, "%lf\n", &beta);
    fscanf(pFile, "%lf\n", &dt);
    fscanf(pFile, "%i\n", &N_samp);
    fscanf(pFile, "%lf\n", &quad);
    fscanf(pFile, "%lf\n", &cub);
    fscanf(pFile, "%lf\n", &quart);
    
    fclose (pFile);
    
    cout << P << endl;
    cout << N << endl;
    cout << corr_length << endl;
    cout << beta << endl;
    cout << dt << endl;
    cout << N_samp << endl;
    cout << quad << endl;
    cout << cub << endl;
    cout << quart << endl;
    
    double** O = new double*[P];
    double** O_T = new double*[P];
    for (int i = 0 ; i < P ; ++i)
    {
        O[i] = new double[P];
        O_T[i] = new double[P];
    }
    // Defining matrix element O_mn
    
    for (int n = 0 ; n < P ; ++n)
    {
        O[0][n] = sqrt(1/static_cast<double>(P));
    }
    
    for (int m = 1 ; m < P/2 ; ++m)
    {
        for (int n = 0 ; n < P ; ++n)
        {
            O[m][n] = sqrt(2.0/static_cast<double>(P)) * cos(2 * M_PI * m * (n+1) / P);
        }
    }
    
    for (int n = 0 ; n < P ; ++n)
    {
        if (n % 2 == 0)
        {
            O[P/2][n] = -sqrt(1/static_cast<double>(P));
        }
        if (n % 2 == 1)
            O[P/2][n] = sqrt(1/static_cast<double>(P));
    }
    
    for (int m = P/2 + 1 ; m < P ; ++m)
    {
        for (int n = 0 ; n < P ; ++n)
        {
            O[m][n] = sqrt(2.0/static_cast<double>(P)) * sin( 2 * M_PI * m * (n+1) / P);
        }
    }
    
    for (int i = 0 ; i < P ; ++i)
    {
        for (int j = 0 ; j < P ; ++j)
        {
            O_T[i][j] = O[j][i];
        }
    }
    
    double* omega = new double[P];
    for (int k = 0 ; k < P ; ++k)
    {
        omega[k] = (2 / (beta_P * h_bar)) * sin(k * M_PI / P);
    }
    
    double* x_normal = new double[P];
    double* x_bead = new double[P];
    double* p_normal = new double[P];
    double* p_bead = new double[P];
    
    double** centroid = new double*[N_samp];
    // Centroid position
    // Store for every run (easier to handle later)
    for (int i = 0 ; i < N_samp ; ++i)
    {
        centroid[i] = new double[N];
    }
    
    double** correlation = new double*[N_samp];
    for (int i = 0 ; i < N_samp ; ++i)
    {
        correlation[i] = new double[corr_length];
    }
    double *corr_avg = new double[corr_length];
    
    // ---------------------------------------------------------------
    for (int rep = 0 ; rep < N_samp ; ++rep)
    {
        // Initialization
        for (int i = 0 ; i < P ; ++i)
        {
            x_normal[i] = 1.0;
            p_normal[i] = momentum_dist();
        }
        
        transform(x_bead, O_T, x_normal);
        transform(p_bead, O_T, p_normal);
        
        centroid[rep][0] = 0.0;
        for (int i = 0 ; i < P ; ++i)
        {
            centroid[rep][0] = centroid[rep][0] + x_bead[i] / P;
        }
        
        // Time propagation
        for (int time = 0 ; time < (N-1) ; ++time)
        {
            // print_centroid(x_bead);
            
            bead_mom_update_half(p_bead, x_bead);
            
            transform(x_normal, O, x_bead);
            transform(p_normal, O, p_bead);
            
            normal_update(p_normal, x_normal, p_normal, x_normal, omega);
            
            transform(x_bead, O_T, x_normal);
            transform(p_bead, O_T, p_normal);
            
            bead_mom_update_half(p_bead, x_bead);
            
            centroid[rep][time+1] = 0.0;
            for (int i = 0 ; i < P ; ++i)
            {
                centroid[rep][time+1] = centroid[rep][time+1] + x_bead[i] / (P);
            }
        }
        
        /*
        for (int i = 0 ; i < N ; ++i)
        {
            cout << centroid[rep][i] << endl;
        }
        cout << endl;
        */
        
        cout << "Run " << rep + 1 << " done out of " << N_samp << endl;
        cout << endl;
    }
    // ---------------------------------------------------------------
    
    /*
     // Exporting the first centroid trajectory
     FILE * pFile;
     pFile = fopen ("centroid.txt", "w");
     for (int i = 0 ; i < N ; ++i)
     {
     fprintf(pFile, "%lf    %lf\n", i*dt, centroid[0][i]);
     }
     fclose (pFile);
     */
    
    // Calculating correlation
    
    cout << "Calculating position ACF..." << endl;
    
    for (int i = 0 ; i < N_samp ; ++i)
    {
        correlation[i] = calc_corr(centroid[i]);
        cout << "Run " << i + 1 << " done out of " << N_samp << endl;
        cout << endl;
    }
    
    double normalizing = 0.0;
    for (int j = 0; j < N_samp ; ++j)
    {
        normalizing = normalizing + correlation[j][0] / N_samp;
    }
    
    cout << "Averaging the position ACF..." << endl;
    for (int i = 0 ; i < corr_length ; ++i)
    {
        corr_avg[i] = 0.0;
        for (int j = 0; j < N_samp ; ++j)
        {
            corr_avg[i] = corr_avg[i] + correlation[j][i] / (N_samp * normalizing);
        }
    }
    
    for (int i = 0 ; i < corr_length ; ++i)
    {
        cout << corr_avg[i] << endl;
    }
    
    cout << "Exporting the position ACF..." << endl;
    FILE * pcorr;
    pcorr = fopen ("corr.txt", "w");
    for (int i = 0 ; i < corr_length ; ++i)
    {
        fprintf(pcorr, "%lf    %lf\n", i*dt, corr_avg[i]);
    }
    fclose (pFile);
     
    delete[] x_normal;
    delete[] x_bead;
    delete[] p_normal;
    delete[] p_bead;
    delete[] omega;
    
    for (int i = 0 ; i < P ; ++i)
    {
        delete[] O[i];
        delete[] O_T[i];
    }
    delete[] O;
    delete[] O_T;
    
    for (int i = 0 ; i < N_samp ; ++i)
    {
        delete[] centroid[i];
        delete[] correlation[i];
    }
    delete[] centroid;
    delete[] correlation;
    
    return 0;
}
