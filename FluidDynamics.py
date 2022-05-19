# -*- coding: utf-8 -*-
"""
Created on Sun May 30 10:22:54 2021

@author: Benjamin Luo
"""
import math

def main():
    
    properties = {"NOX":[787.7069930752898, 1.47*10**(-5)],"GOX":[0, 2.04*10**(-5)],"LOX":[]}
    
    roughness = {"Stainless Steel":3.5*10**(-6),"Aluminium":1.5*10**(-6),"Copper":1.5*10**(-6)}
    
    
    print("This program assumes a circular cross section pipe and all metric IS units")
    print("Temperature is assumed to be 293K")
    with open("fluidData.txt","r") as f:
        data = f.readlines()
    length = float(data[0].split(": ")[-1])
    pipe_dia = float(data[1].split(": ")[-1])/1000
    pipe_material = data[2].split(": ")[-1].strip("\n")
    flow_rate = float(data[3].split(": ")[-1])
    fluid = data[4].split(": ")[-1].strip("\n")
    fluid_density = properties[fluid][0]
    fluid_dviscosity = properties[fluid][1]
    flow_velocity = flow_rate/(fluid_density*math.pi*(pipe_dia/2)**2)
    Re = fluid_density*flow_velocity*pipe_dia/fluid_dviscosity
    k = roughness[pipe_material]
    
    if Re < 2300:
        print("Flow is laminar")
        friction_c = 64*fluid_dviscosity/(pipe_dia*fluid_density*flow_velocity)
    else:
        print("Flow is turbulent")
        friction_c = guess_lambda(1/(2*math.log10(k/(3.72*pipe_dia)))**2,Re,k,pipe_dia)
        # friction_c = 1/(2*math.log10(k/(3.72*pipe_dia)))**2 #Approximation for large Re
    
    pressure_loss = friction_c*length/pipe_dia*fluid_density*flow_velocity**2/2
    print("Reynold's number: ", Re)
    print("Fluid density:",fluid_density, "kg m^-3")
    print("Friction coefficient:", friction_c)
    print("Flow velocity:",flow_velocity,"m s^-1")
    print("Pressure Loss:", pressure_loss, "Pa,",pressure_loss/10**5,"bar")
            
def guess_lambda(guess_c,Re,k,pipe_dia):
    guess = -2*math.log10(2.51/(Re*guess_c**0.5)+k/(3.72*pipe_dia))-guess_c**(-0.5)
    # print("Guess:", guess)
    # print("Friction c.:", guess_c)
    if guess > 0.0002:
        return guess_lambda(guess_c-0.000001,Re,k,pipe_dia)
    elif guess <-0.0002:
        return guess_lambda(guess_c+0.000001,Re,k,pipe_dia)
    else:
        return round(guess_c, 6)
            
    

main()