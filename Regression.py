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
a = 0.236 #0.17
n = 0.5
length = 0.1 # 0.7
fuel_density = 1010 #900
throat_A =  0.000613598 #0.0000224 # 
v_e = 1965 #2338.664

def g_ox(ox_flowRate):
    return ox_flowRate*10**3/(math.pi*portR**2)

def fuelRegression(g_ox):
    return a*g_ox**n/10**3

def fuel_flowRate(fuelRegression):
    return fuelRegression*2*math.pi*portR/100*length*fuel_density

def getPC(ox_flowRate):
    mass_flow = fuel_flowRate(fuelRegression(g_ox(ox_flowRate)))+ox_flowRate
    p_c = 0.5*mass_flow*math.sqrt(1.15*8.3144621*3330.72/0.022897959)/(throat_A*1.15*math.sqrt((2/2.15)**(2.15/0.15)))/10**5
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