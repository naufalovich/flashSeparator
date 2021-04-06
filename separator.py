"""
Created on Sun Apr 04 2021 17:00:15 

@author: naufalovich
"""
# This program perform calculating Bubble/ Dew point of binary mixture
# using phasepy libray and peng-robinson equation of state

# import necessary library
import numpy as np
from phasepy import component, mixture, preos
from phasepy.equilibrium import flash
import scipy.linalg

# given data
F = 1000 # molar flow. kmol/h
T = 85. + 273.15  # vessel temperature, converted to K
P = 1.01 # vessel pressure, bar
Z = np.array([0.75, 0.25]) # overall molar fraction component

# define component thermodynamical properties
benzene = component(name='benzene', Tc=562.2, Pc=48.98, Zc=0.271, Vc=259.0, w=0.210,
                        Ant=[13.7819, 2726.81, 217.572],
                        GC={'CH=C':6})

toluene = component(name='toluene', Tc=591.8, Pc=41.06, Zc=0.264, Vc=316, w=0.262,
                      Ant=[13.9320, 3056.96, -217.625],
                      GC={'CH=C':6, 'CH3':1})

# setting up eos
mix = mixture(benzene, toluene) #mix given component
mix.unifac() # using dortmund modified unifac mixing rule
eos = preos(mix, 'mhv_unifac') # peng robinson equation of state

# initial guess
x0 = np.array([0.4, 0.6]) # liquid molar fraction to guess
y0 = np.array([0.2, 0.8]) # vapor molar fraction to guess

# start calculation
sep = flash(x0, y0, 'VL', Z, T, P, eos) # phase compositions, vapor phase fraction
y1 = (sep[0])[0] # benzene vapor fraction, indexing from tuple
y2 = (sep[0])[1] # toluene vapor fraction
x1 = (sep[1])[0] # benzene liquid fraction
x2 = (sep[1])[1] # benzene liquid fraction

# compute vapor flow and liquid flow using matrix
A = np.array ([[y1, x1],[y2, x2]])
B = np.array ([F*Z[0], F*Z[1]])

# now solve the matrix
C = scipy.linalg.solve(A,B) # solver for linear algebra

# print all solution
print ("-----------------------------------------------")
print ("Vapor Flow (V) :", C[0], "kmol/h")
print ("     benzene fraction in vapor (y1) :", y1)
print ("     toluene fraction in vapor (y2) :", y2)
print ("-----------------------------------------------")
print ("Liquid Flow (L) :", C[1], "kmol/h")
print ("     benzene fraction in liquid (x1) :", x1)
print ("     toluene fraction in liquid (x2) :", x2)
