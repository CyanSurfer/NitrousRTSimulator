# -*- coding: utf-8 -*-
"""
Created on Fri May 14 09:46:31 2021

@author: Benjamin Luo
"""

import math

#Initialise Variables
#Critical Pressure = 72.51 bar
#Critical density = 452 kg m^-3
#g_ox in g cm^-2, fuelRegression in ms^-1
portR = 1 #2.5
a = 0.0766 #0.17
n = 0.795
length = 0.1 # 0.7
fuel_density = 1010 #900
throat_A =  0.00012327 #0.000095 # #0.000613598 #0.0000224 # 
v_e = 1965 #2338.664
gamma = 1.2018
combustion_T = 2886.92
s_gas_constant = 300.86

def g_ox(ox_flowRate):
    return ox_flowRate*10**3/(math.pi*portR**2)

def fuelRegression(g_ox):
    return a*g_ox**n/10**3

def fuel_flowRate(fuelRegression):
    return fuelRegression*2*math.pi*portR/100*length*fuel_density

def getPC(ox_flowRate):
    mass_flow = fuel_flowRate(fuelRegression(g_ox(ox_flowRate)))+ox_flowRate
    p_c = mass_flow*math.sqrt(combustion_T*s_gas_constant/gamma)/(throat_A*gamma*math.sqrt((2/(gamma+1))**((gamma+1)/(gamma-1))))/10**5
    print("Combustion Chamber Pressure / bar:", p_c)
    return p_c

def regressFuel(ox_flowRate):
    global portR
    portR+=fuelRegression(g_ox(ox_flowRate))*100
    
def getThrust(oxi_flowRate):
    return v_e*(fuel_flowRate(fuelRegression(g_ox(oxi_flowRate)))+oxi_flowRate)

def update_ve(pc):
    global v_e

    new_ve = ...

    v_e = new_ve