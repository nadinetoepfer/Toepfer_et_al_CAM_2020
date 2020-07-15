#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt

import helping_functions as hf


def plot_generated_data(parameter_file, th_dict, sigma,constants):
    
    """
    Initiated by main()
    
    parameters:
    ------
     - parameter file
     - th_dict
     - sigma (parameter for gauss curve that makes the curve match a certain day length)
   
    """

    # light curve
    min_light = 0
    mid_day = float(constants['mid_day'])
    Delta_light = float(th_dict['Max_Light'])
    shift_light = 0
    
    # temperature curve
    minT = float(th_dict['Min_T'])
    maxT = float(th_dict['Max_T'])

    # humidity curve
    max_humidity = float(th_dict['Max_H'])
    min_humidity = float(th_dict['Min_H'])

    # calculations
   
    mean_T = (maxT + minT)/2
    A_T = (maxT - minT)/2
    d_T = float(th_dict['skew_T'])
    per_T = float(th_dict['period_T'])
    
    mean_H = (max_humidity + min_humidity)/2  
    A_H =  - (max_humidity - min_humidity)/2 # try reversing amplitude to get reverse curve to T
    d_H = float(th_dict['skew_H'])
    per_H = float(th_dict['period_H'])
   
    deltaT  = maxT-minT
    
    xx=np.linspace(0, 23, 100)
    
    g_humidity=list()
    for i in xx:
        g_humidity.append(hf.skewed_sinus(i, A_H, per_H, mean_H, d_H))
    
    g_temperature=list()
    for i in xx:
        g_temperature.append(hf.skewed_sinus(i, A_T, per_T, mean_T, d_T))
    
    g_light=list()        
    for i in xx:
        g_light.append(hf.gaussian(i,min_light, Delta_light, mid_day, sigma, shift_light))
      
    
    fig, ax1 = plt.subplots()
    
    ln1=ax1.plot(xx, g_light,'y',label='Light',linewidth=2.8)
    ax1.set_ylim([0,1.1* Delta_light])
    ax1.set_ylabel('Light intensiy in $\mu mol/m^2 \cdot s$', fontsize=14)
    lns = ln1
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0,frameon=False, fontsize=13)
    fig.tight_layout()
    plt.show()
    fig.savefig(os.path.join("Results/Light.png"),dpi=600)

    fig, ax1 = plt.subplots(figsize=(5,2)) 
    
    ax2 = ax1.twinx()
    ax1.set_ylim([minT-2,minT+deltaT+2])

    ax2.set_ylim([min_humidity-0.1,1.1])
    plt.rcParams['xtick.labelsize']=12
    plt.rcParams['ytick.labelsize']=12

    
    ln1=ax1.plot(xx, g_temperature,'orangered',label='Temperature',linewidth=2.8)
    ln2=ax2.plot(xx,g_humidity,'steelblue',label='Relative Humidity',linewidth=2.8)
    
    lns = ln1+ln2
    labs = [l.get_label() for l in lns]
    fig.tight_layout()
    plt.show()
    fig.savefig(os.path.join("Results/Temp_Hum.png"),dpi=600)

    
