# -*- coding: utf-8 -*-
"""
Created on Thu May 13 17:43:10 2021

@author: Benjamin Luo
"""

import pylab

def densityV(temp):
    temp/=309.57
    term = 1/temp-1
    return 452*pylab.exp(-1.009*term**(1/3)-6.28792*term**(2/3)+7.50332*term-7.90463*term**(4/3)+0.629427*term**(5/3))

def densityL(temp):
    temp/=309.57
    term = 1-temp
    return 452*pylab.exp(1.72328*term**(1/3)-0.8395*term**(2/3)+0.5106*term-0.10412*term**(4/3))

def specificHeatL(temp):
    temp/=309.57
    term = 1-temp
    return 2.49973*10**3*(1+0.023454*term**(-1)-3.80136*term+13.0945*term**2-14.5180*term**3)

def latentHeatV(temp):
    temp/=309.57
    term = 1-temp
    h_l = -200 + 116.043*term**(1/3)-917.225*term**(2/3)+794.779*term-589.587*term**(4/3)
    h_g = -200 + 440.055*term**(1/3)-459.701*term**(2/3)+434.081*term-485.338*term**(4/3)
    return (h_g-h_l)*10**3

def vapourPressure(temp):
    temp/=309.57
    term = 1-temp
    return 72.51*pylab.exp((-6.71893*term+1.35966*term**(3/2)-1.3779*term**(5/2)-4.051*term**5)/temp)

# print(densityL(303))