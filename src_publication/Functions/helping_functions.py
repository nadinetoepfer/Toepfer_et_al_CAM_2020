#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import math

def gaussian(x, level, delta, mean, sigma, shift):
    
    """
    - generates normal distributed curves to simulate light
    - curves are shifted so they form a continious curve over the 24h simulation periode
    
    parameters:
    ---------
    - level: baseline value for the curve
    - delta: hights of the curve
    - mean: mid-day
    - shit of curve from mid day towards afternoon
    - sigma: widths of the cuve - calculated to match sunrise 
    
    """

    if x > (24 - (mean - shift)):
        value = level + delta * np.exp(-np.power(x - 24 - (mean + shift), 2.) / (2 * np.power(sigma, 2.)))
    else:
        value = level + delta * np.exp(-np.power(x - (mean + shift), 2.) / (2 * np.power(sigma, 2.)))
    return value


def skewed_sinus(x, A, period, mean, d, f = (2 * math.pi)/24):
    
    """
    - generates a skewed sinus curves to simulate  temperature
    
    parameters:
    ---------
    - A: Amplitude = (Tmax - Tmin)/2 # determines the deviation from mean
    - m: mean temp = (Tmax + Tmin)/2  # sets the temperature around which the sin curve is zentered
    - f: frequency = (2 * math.pi)/24 # periode i.e. 24 hours
    - p : phase
    - d: skew parameter
    - d and p both influence the shape of the curve in interdependency and need to be adjusted together
    
    """
    
    value = mean + A * math.sin(f * ((x - period) + d * math.sin(f *(x - period))/2))
    return value

