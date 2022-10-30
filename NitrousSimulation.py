# -*- coding: utf-8 -*-
"""
Created on Thu May 13 16:23:30.8 2021

@author: Benjamin Luo
"""
import pylab
import math
import NitrousProperties as NP
import Regression as RG

#variables initialisation
#------------------------
#masses in kg
#densities in kg m**-3
#specific heat capacity in J kg**-1 K**-1
#temperature in Kelvin
#volume in m**3
#Pressure in bars
#time in seconds
#Area in m**2
#Specific gas constant in J kg^-1 K^-1
#Critical Z = 0.28
#Injector diameter in mm


mass_v = 0
mass_l = 0
mass_vapourised = 0
density_v = 0
density_l = 0
c_liquid = 0
temp = 300
V_tank = 0.0034
oxi_fill = 0.85

P_tank = 0
P_manifold = 0
P_injector = 0
P_loss = 0

Z = 0
Z_2 = 0
gamma = 1.3
R=188.91

time = 0
time_step = 0.1

K = 0.8 #1.0
A_injector = 0
injector_diameter = 1.75
N = 1
D_loss = 0

burnout = False

initial_vTemp = 0
initial_vDensity = 0
initial_vPressure = 0
initial_vMass = 0
thrust = 0

rocketMass = 30
g = 9.81

Error = False

def initialise():
    
    global density_v,density_l,mass_v,mass_l,A_injector,D_loss,P_tank,V_tank, temp, Error
    try:

        temp = float(input("Enter initial Nitrous temperature in Kelvins:"))

        density_v = NP.densityV(temp)
        density_l = NP.densityL(temp)

        V_tank = float(input("Enter the Nitrous Tank Volume in litres:"))

        if V_tank <= 0:
            total_mass = float(input("Enter the Mass of Liquid Nitrous Oxide needed in kg:"))
            
            mass_l = total_mass/(1+density_v/density_l*(1/oxi_fill-1))
            mass_v = total_mass - mass_l
            print(mass_l, mass_v)
            V_tank = mass_l/density_l + mass_v/density_v
            
            print("Total Volume with "+str(1-oxi_fill)+" ullage", V_tank*1000, "litres")
            
        else:

            V_tank = V_tank/1000
            
            mass_v = (1-oxi_fill)*V_tank*density_v
            mass_l = oxi_fill*V_tank*density_l
            
            print("Total mass:", mass_v+mass_l)


        A_injector = math.pi*((injector_diameter/1000)**2)/4
        D_loss=K/(N*A_injector)**2
        #print("hello",D_loss)
        P_tank = NP.vapourPressure(temp)
        # print(mass_v,mass_l)
        
    except:
        
        print("There was an error in the input.")
        
        Error = True

def injector_flowRate(C_d, A, density, P_manifold, P_injector):
    
    return C_d*A*math.sqrt(2*density*(P_manifold-P_injector)*10**5)
    

def guessZ(Z_guess):
    T_2 = ((Z_guess*mass_v)/(Z*initial_vMass))**(gamma-1)*initial_vTemp
    # print("Initial pressure",initial_vPressure)
    # print("new Temp",T_2)
    P_2 = (T_2/initial_vTemp)**(gamma/(gamma-1))*initial_vPressure
    # print("new pressure",P_2)
    # Z_prime = (P_2/initial_vPressure)**(1/gamma)*Z*initial_vMass/mass_v
    Z_prime = 1-24/2417*P_2*1.01325
    # print("Iterated Z",Z_prime,Z_guess,T_2,P_2)
    if abs(Z_prime-Z_guess)/Z_prime>0.01:
        return guessZ((Z_guess+Z_prime)/2)
    else:
        return Z_guess,T_2,P_2

def loop():

    global density_v,density_l,mass_v,mass_l,A_injector,D_loss,time, \
        time_step,temp,P_tank,P_manifold,P_injector,burnout,mass_vapourised, thrust, rocketMass, \
        Error
        
    P_manifold = P_tank-P_loss
    heatR = mass_vapourised*NP.latentHeatV(temp)
    deltaT = -heatR/(mass_l*NP.specificHeatL(temp))
    # print("Temperature Change",deltaT)
    temp += deltaT
    density_v = NP.densityV(temp)
    density_l = NP.densityL(temp)
    #ml_dot = (2*density_l*(P_manifold-P_injector)*10**5/D_loss)**0.5
    ml_dot = injector_flowRate(K, A_injector, density_l, P_manifold, P_injector)
    rocketMass-=(ml_dot+RG.fuel_flowRate(RG.fuelRegression(RG.g_ox(ml_dot))))*time_step
    thrust = RG.getThrust(ml_dot)
    RG.regressFuel(ml_dot)
    P_injector = RG.getPC(ml_dot)
    if P_injector > P_manifold:
        print("Combustion chamber pressure:", P_injector)
        print("Combustion chamber pressure too high! Back flow expected")
        Error = True
        return None
    # print("rate of nitrous discharge",ml_dot)
    ml_old = mass_l-ml_dot*time_step
    ml_new = (V_tank-(ml_old+mass_v)/density_v)/(1/density_l-1/density_v)
    mass_vapourised = ml_old-ml_new
    print("Mass Vapourised", mass_vapourised)
    # print("mass flow rate:",ml_dot)
    # print("Chamber Pressure",P_injector)
    print(mass_v/density_v)
    if mass_vapourised<0 or mass_l<0:
        # ml_old becomes negative
        
        print("Burnout!")
        burnout = True
        # mass_v += mass_l
        
        # print("Temperature Change",deltaT)
        # flashVapourise()
    else:
        mass_l = mass_l - mass_vapourised - ml_dot*time_step
        mass_v += mass_vapourised
        time+=time_step
        
    print("Mass of liquid Nitrous",mass_l)
    print("Temperature",temp)   
    
    
    
    P_tank = NP.vapourPressure(temp)
    #print("REgression Rate:", RG.fuelRegression(RG.g_ox(ml_dot)),RG.portR)
    # print(temp)

# def flashVapourise():
#     global temp
    
#     for i in range(int(math.ceil(mass_l/0.05))):
#         heatR = min(mass_l-i*0.05,0.05)*NP.latentHeatV(temp)
#         print(heatR)
#         deltaT = -heatR/((mass_l-(i+1)*0.05)*NP.specificHeatL(temp))
#         temp+=deltaT

def vapourPhase():
    
    global density_v,density_l,mass_v,mass_l,A_injector,D_loss,time, \
        time_step,temp,P_tank,P_manifold,P_injector,burnout,mass_vapourised,Z, \
            initial_vTemp,initial_vDensity,initial_vPressure,initial_vMass,Z_2, thrust, rocketMass
            

    
    Z = 1-24/2417*P_tank*1.01325
    
    P_manifold=P_tank-P_loss
    print("mass of liquid",mass_l)
    print("mass of vapour",mass_v)

    density_v = mass_v/V_tank
    initial_vDensity = density_v
    initial_vPressure = P_tank#Z*density_v*temp*R/10**5
    # P_tank = initial_vPressure
    print("Pressures:",P_tank,initial_vPressure)
    initial_vTemp = temp
    initial_vMass = mass_v
    
    
    #mv_dot = (2*density_v*(P_manifold-P_injector)*10**5/D_loss)**0.5
    mv_dot = injector_flowRate(K, A_injector, density_v, P_manifold, P_injector)
    rocketMass-= (mv_dot+RG.fuel_flowRate(RG.fuelRegression(RG.g_ox(mv_dot))))*time_step
    thrust = RG.getThrust(mv_dot)
    RG.regressFuel(mv_dot)
    P_injector = RG.getPC(mv_dot)
    mass_v -= mv_dot*time_step
    print("first mass drop:",mass_v)
    Z_2,temp,P_tank = guessZ(1-24/2417*P_tank*1.01325)
    density_v = mass_v/V_tank
    time+=time_step
    
def vapourLoop():
    global P_tank,P_manifold,P_injector,density_v,temp,mass_v,time,Z_2, thrust, rocketMass
    
    
    P_manifold=P_tank-P_loss
    # print("density",density_v,initial_vPressure,P_manifold, P_injector)
    #mv_dot = (2*density_v*(P_manifold-P_injector)*10**5/D_loss)**0.5
    mv_dot = injector_flowRate(K, A_injector, density_v, P_manifold, P_injector)
    rocketMass-= (mv_dot+RG.fuel_flowRate(RG.fuelRegression(RG.g_ox(mv_dot))))*time_step
    thrust = RG.getThrust(mv_dot)
    RG.regressFuel(mv_dot)
    print("Vapour mass flow:",mv_dot)
    Z_2,temp,P_tank = guessZ(1-24/2417*P_tank*1.01325)
    density_v = mass_v/V_tank
    P_injector = RG.getPC(mv_dot)
    P_tank = Z_2*density_v*R*temp/10**5
    mass_v -= mv_dot*time_step
    
    time+=time_step

def plot(p_tank, thrustPlot, rocket_mass):
    
    rocket_vel = []
    altitude = []
    acceleration = []
    totalImpulse = 0
    
    xVals = pylab.array([x*time_step for x in range(0,round(time/time_step)+1)])
    rocket_vel = [0.0 for x in range(0,round(time/time_step)+1)]
    # print(len(xVals))
    p_tank = pylab.array(p_tank)
    
    
    # thrustPlot = pylab.array(thrustPlot)
    thrustPlot = pylab.array(thrustPlot)
    for i in range(len(thrustPlot)):
        totalImpulse+=thrustPlot[i]*time_step
        acceleration.append(thrustPlot[i]/rocket_mass[i])
        if i>0: 
            rocket_vel[i]+=0.5*(acceleration[-1]+acceleration[-2])*time_step+rocket_vel[i-1]
            #print(rocket_vel[i],acceleration[-1])
            altitude.append(0.5*(rocket_vel[i]+rocket_vel[i-1])*time_step+altitude[-1])
        else:
            altitude.append(0.0)
    
    #coasting
    while rocket_vel[-1]>0:
        acceleration.append(-g)
        rocket_vel.append(rocket_vel[-1]-g*time_step)
        altitude.append(0.5*(rocket_vel[-1]+rocket_vel[-2])*time_step+altitude[-1])
    
    rocket_vel = pylab.array(rocket_vel)
    acceleration = pylab.array(acceleration)
    
    altitudeTime = pylab.array([x*time_step for x in range(len(altitude))])
    # print(rocket_vel)
    # pressure_plot, thrust_plot, altitude_plot = pylab.subplots(1,3)
    
    fig, ax = pylab.subplots()

    
    # pylab.plot(xVals,yVals,"-k",label="Pressure")
    ax.plot(altitudeTime, acceleration,"r.",label="Acceleration")
    ax.set_xlabel("Time/s")
    ax.set_ylabel("Acceleration/ms^-2")
    ax2=ax.twinx()
    ax2.plot(altitudeTime, pylab.array(altitude), "-b",label="altitude")
    ax2.set_ylabel("Altitude/m")
    fig.show()
    
    pylab.figure()
    pylab.plot(xVals, p_tank,"r+",label="RT_P")
    pylab.xlabel("Time/s")
    pylab.ylabel("Run Tank Pressure / bar")
    pylab.show()
    
    print("Average Thrust:", sum(thrustPlot)/len(thrustPlot))
    print("Altitude:", altitude[-1])
    print("Total Impulse:", totalImpulse)
    print("Final Rocket Mass:", rocket_mass[-1])
    # pressure_plot.plot(xVals,yVals,"k+",label="Pressure")
    # thrust_plot.plot(xVals, thrustPlot,"r+",label="Thrust")
    # altitude_plot.plot(altitudeTime, pylab.array(altitude), "b+",label="altitude")
    
    # pressure_plot.xlabel("Time/s")
    # pressure_plot.ylabel("Pressure/bar")
    # thrust_plot.xlabel("Time/s")
    # thrust_plot.ylabel("Thrust/N")
    # altitude_plot.xlabel("Time/s")
    # altitude_plot.ylabel("Altitude/m")
    
    fig.legend()

def main():
    
    p_tank = []
    thrustPlot = []
    temperature = []
    rocket_mass = []
    
    initialise()
    
    if Error:
        return None
    
    while not burnout:
        loop()
        if not Error:
            temperature.append(temp)
            p_tank.append(P_tank)
            thrustPlot.append(thrust)
            rocket_mass.append(rocketMass)
        else:
            break
    
    if not Error:
        vapourPhase()
        print("Initial vapour Mass",mass_v)
        p_tank.append(P_tank)
        temperature.append(temp)
        thrustPlot.append(thrust)
        rocket_mass.append(rocketMass)
        
        #print(P_tank,Z,Z_2,temp,mass_v)
        #print(len(yVals))
        while mass_v>0.01:
            vapourLoop()
            p_tank.append(P_tank)
            temperature.append(temp)
            thrustPlot.append(thrust)
            rocket_mass.append(rocketMass)
            print(RG.portR)
        
        plot(p_tank, thrustPlot, rocket_mass)
    
    
if __name__ == "__main__":    
    main()
