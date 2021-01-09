#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 15:05:10 2021

Helper functions for main

@author: ebalestr
"""
import numpy as np




def a_prime(HR, FA):
    
    return .5 + np.sign(HR-FA)*(((HR-FA)**2 + np.abs(HR-FA)) / 
                                  (4*np.max([HR, FA]) - 4*HR*FA))



