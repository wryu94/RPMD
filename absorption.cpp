//
//  main.cpp
//  RPMD_absorption_spectrum_calculation
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

int main()
{
    // Can calculate the spectrum as a part of RPMD, but too lazy & already started writing this, so...
    int N = 10000; // How long the RPMD position ACF are
    int P = N; // How long the cosine transform are
    
    double *correlation = new double[N];
    double *time = new double[N];
    
    cout << "Reading in position ACF..." << endl;
    cout << endl;
    
    FILE *pFile;
    pFile = fopen ("position_ACF.txt", "r"); // Reading in the correlation function
    for (int i = 0 ; i < N ; ++i)
    {
        fscanf(pFile, "%lf", &time[i]);
        fscanf(pFile, "%lf", &correlation[i]);
    }

    double dw = (2 * M_PI) / time[N-1]; // Smallest unit of frequency
    double dt = time[1];
    
    cout << "Calculating the absorption spectrum..." << endl;
    cout << endl;
    
    double *absorbance = new double[P];
    for (int k = 0 ; k < P ; ++k)
    {
        absorbance[k] = 0.0;
        for (int time = 0 ; time < N-1 ; ++time)
        {
            absorbance[k] = absorbance[k] + (k * k * dw * dw) * (dt / 2) * (correlation[time] * cos((k * dw) * (time * dt)) + correlation[time+1] * cos((k * dw) * ((time+1) * dt)));
        }
        cout << static_cast<double>(k+1)*100/static_cast<double>(P) << "% done..." << endl;
    }
    
    cout << "Exporting the absorbance spectrum..." << endl;
    cout << endl;
    
    ofstream myfile;
    myfile.open ("absorbance.txt");
    
    for (int k =0 ; k < P ; ++k)
    {
        myfile << k * dw << " " << absorbance[k] << endl;
    }
    
    myfile.close();
    
    delete[] correlation;
    delete[] time;
    delete[] absorbance;
    
    return 0;
}
