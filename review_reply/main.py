#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 15:04:36 2021

@author: ebalestr

The script aims to reply to reviewers for the manuscript to EJN
Additional analyses not present here but provided to reviewers have been performed 
by changeing starting parameters in *.m functions and scripts written for the 
original data analysis and stored in ../main_data_analysis

"""

import numpy as np
from scipy.io import loadmat
from scipy.stats import norm
import matplotlib.pyplot as plt
from helper import a_prime


# load data from mat file
dat = loadmat('../RAW_data_collapsed')['mat_data']
nsubjs = dat.shape[2];

# rows reminder:
# 1 row: load condition
# 2 row: mem equality condition
# 3 row: subj mem response
# 4 row: deltaT
# 5 row: timestamp delta
# 6 row: flash_presence
# 7 row: subj flash response

# in no particular order, response to single reviewers points follows

#%% VWM performance (R1, points 7, 8, 10)

# accuracy (7) , HR and FA of VWM task (8)
isequal_VWM = dat[1, :, :]
reported_equal_VWM = dat[2, :, :]
flashpresent = dat[5, :, :]
SOA_flash = dat[3, :, :]

SOAs_vals = np.unique(SOA_flash)


# define here load masks
load_masks =  {}
for iload in [0, 2, 4]:
    
    load_masks['load_' + str(iload)] = dat[0, :, :] == iload

# for each part subselect load 2 and 4, the current SOA 
# (in which a flash had been presented) and compute metrics

acc_L2 = []
acc_L4 = []
L2_HR = []
L4_HR = []
L2_FA = []
L4_FA = []

VWM_trends_SOA = np.zeros((len(SOAs_vals), nsubjs, 2))

for isubj in range(nsubjs):
    
    this_isequal_VWM = isequal_VWM[:,isubj]
    this_repequal_VWM = reported_equal_VWM[:,isubj]
    L2_mask = load_masks['load_2'][:, isubj]
    L4_mask = load_masks['load_4'][:, isubj]
    flash_present_mask = flashpresent[:, isubj]==1
    
    # accuracy
    match = (this_isequal_VWM == this_repequal_VWM)*1
    acc_L2.append(match[L2_mask].mean())
    acc_L4.append(match[L4_mask].mean())
    
    # SDT    
    L2_HR.append(this_repequal_VWM[(this_isequal_VWM==1) &
                                   L2_mask].mean())
    L2_FA.append(this_repequal_VWM[(this_isequal_VWM==0) &
                              L2_mask].mean())
    L4_HR.append(this_repequal_VWM[(this_isequal_VWM==1) &
                              L4_mask].mean())
    L4_FA.append(this_repequal_VWM[(this_isequal_VWM==0) &
                              L4_mask].mean())
    
    acc_SOA = 0
    for isoa in SOAs_vals:
        
        SOA_mask = SOA_flash[:, isubj] == isoa
        L2_soa_HR = this_repequal_VWM[(this_isequal_VWM==1) &
                                      L2_mask &
                                      SOA_mask &
                                      flash_present_mask].mean()
    
        L2_soa_FA = this_repequal_VWM[(this_isequal_VWM==0) &
                                      L2_mask &
                                      SOA_mask &
                                      flash_present_mask].mean()

        L4_soa_HR = this_repequal_VWM[(this_isequal_VWM==1) &
                                      L4_mask &
                                      SOA_mask &
                                      flash_present_mask].mean()

        L4_soa_FA = this_repequal_VWM[(this_isequal_VWM==0) &
                                      L4_mask &
                                      SOA_mask &
                                      flash_present_mask].mean()

        VWM_trends_SOA[acc_SOA, isubj, 0] = a_prime(L2_soa_HR, L2_soa_FA)
        VWM_trends_SOA[acc_SOA, isubj, 1] = a_prime(L4_soa_HR, L4_soa_FA)

        acc_SOA += 1


# produce plots
plt.figure()

# HR and FAs
ax = plt.subplot(1, 4, 4)

x = np.array([1, 2])

HR_VWM = [np.mean(L2_HR), np.mean(L4_HR)]
FA_VWM = [np.mean(L2_FA), np.mean(L4_FA)]
width = .35

HR_err = [np.std(L2_HR)/np.sqrt(nsubjs), np.std(L4_HR)/np.sqrt(nsubjs)]
FA_err = [np.std(L2_FA)/np.sqrt(nsubjs), np.std(L4_FA)/np.sqrt(nsubjs)]

        
ax.bar(x - width/2, HR_VWM, width, label='HR', yerr=HR_err)
ax.bar(x + width/2, FA_VWM, width, label='FA', yerr=FA_err)

ax.set_xticks(x)
ax.set_xticklabels(['load 2', 'load4'])
ax.legend()
plt.title('VWM HR and FA')

# A prime trend over SOAs
avg_L2 = VWM_trends_SOA[:, :, 0].mean(axis=1)
avg_L4 = VWM_trends_SOA[:, :, 1].mean(axis=1)

err_L2 = VWM_trends_SOA[:, :, 0].std(axis=1)/np.sqrt(nsubjs)
err_L4 = VWM_trends_SOA[:, :, 1].std(axis=1)/np.sqrt(nsubjs)


c2 = np.array([127, 0, 255])/255
c4 = np.array([255, 51, 153])/255


plt.subplot(1, 4, (1,3))
plt.plot(SOAs_vals*10, avg_L2, c=c2, label='load 2')
plt.fill_between(SOAs_vals*10, avg_L2 - err_L2, avg_L2 + err_L2, color=c2, alpha=.15)

plt.plot(SOAs_vals*10, avg_L4, c=c4, label='load 4')
plt.fill_between(SOAs_vals*10, avg_L4 - err_L4, avg_L4 + err_L4, color=c4, alpha=.15)

plt.xlabel('SOA')
plt.ylabel('A prime')
plt.title('VWM sensitivity as a function of flash appearance SOA')
plt.legend()

















































    
    