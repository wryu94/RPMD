# RPMD

There are two source codes contained in this folder: reduced_unit_RPMD.cpp and absorption.cpp. 
The first code runs the RPMD trajectory of a ring polymer moving along 1D potential and calculate the position autocorrelation function (ACF) from said trajectories. The second code uses the position ACF from the first code to calculate the absorbance spectrum. 

The theoretical backgroun on RPMD simulation can be inferred from the original RPMD paper from 2004 by Manolopoulos. The gist is that you set the ficticious mass equal to the physical mass such that motion of all the beads become all physical. The time propagation step in the RPMD code is from the 2010 PILE paper also from Manolopoulos, which uses the classical Liouvillian operator formulation. 

Few things to note about the RPMD code is that you're sampling momentum of the normal modes using the elevated temperature beta/P, not regular beta. Another thing to be careful about is the unit implemented in the code. There are no explicit units for the parameters, so if one wants to translate this into simulation of more realistic system, one would need to think more specifically about converting the units, etc.

One can obtain absorption spectrum by taking the cosine transform of the position autocorrelation function. The relationship can be written as: Abs(w) = w^2 * integrate(cos(w*t) <x(0)x(t)>,t,0,infinity). The code performs the integral by simple trapezoid rule. 
