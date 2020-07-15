#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd
import math
import seaborn as sns

# define the nesesary parameters
def roundup(x):
    return int(math.ceil(x / 10.0)) * 10

def rounddown(x):
    return int(math.floor(x / 10.0)) * 10

'''
get result files
'''

path = '<add path to results file>'

file_ICDH = '<add file name here>'
file_noICDH = '<add file name here>'
    
file_name_ICDH = path + file_ICDH
file_name_noICDH =  path + file_noICDH

df_ICDH = pd.read_csv(file_name_ICDH, header=0, sep=',')
df_noICDH = pd.read_csv(file_name_noICDH, header=0, sep=',')

# turn values to integer for nicer plot axis
df_ICDH['MaxT'] = df_ICDH['MaxT'].round(0).astype(int)
df_ICDH['MinT'] = df_ICDH['MinT'].round(0).astype(int)

df_noICDH['MaxT'] = df_noICDH['MaxT'].round(0).astype(int)
df_noICDH['MinT'] = df_noICDH['MinT'].round(0).astype(int)

# which pareto step to plot
pareto_08 = 'Water loss_0.8'
pareto_1 = 'Water loss_1.0'


'''
organize data
'''

df_p_ICDH = df_ICDH
df_p_noICDH = df_noICDH

# sort by MinH and MaxT
M_1_noICDH=df_p_noICDH[['MinHumidity', 'MaxT', pareto_1]] # noICDH at pareto_1 
M_1_noICDH_abs = M_1_noICDH.pivot(index = 'MinHumidity', columns ='MaxT',values = pareto_1)

M_08_ICDH=df_p_ICDH[['MinHumidity', 'MaxT', pareto_08]] # ICDH at Pareto 0.8
M_08_ICDH_abs = M_08_ICDH.pivot(index = 'MinHumidity', columns ='MaxT',values = pareto_08)

M_08_noICDH=df_p_noICDH[['MinHumidity', 'MaxT', pareto_08]] # noICDH at pareto 0.8
M_08_noICDH_abs = M_08_noICDH.pivot(index = 'MinHumidity', columns ='MaxT',values = pareto_08)

# some scaling factors, can be adjusted
a = max(M_1_noICDH_abs/1000)
b = 5


'''
caluculate values to plot
'''

# absolute water loss C3 scenario of noICDH (because this might have slightly higher water loss) - CAM with ICDH at 80% 
Diff_abs = M_1_noICDH_abs.subtract(M_08_ICDH_abs, fill_value=0) 

# relative water loss C3 scenario of noICDH - CAM with ICDH at 80% 
Diff_rel = M_1_noICDH_abs.subtract(M_08_ICDH_abs, fill_value=0) *100/M_1_noICDH_abs

# absolute contribution of ICDH with respect to noICDH at 80 %
Diff_abs_contrib = M_08_noICDH_abs.subtract(M_08_ICDH_abs, fill_value=0) 

# realative contribution of ICDH with respect to noICDH at 80 % from C3 (100%) noICDH 
Diff_rel_contrib = M_08_noICDH_abs.subtract(M_08_ICDH_abs, fill_value=0)*100/M_1_noICDH_abs #  (noI_80 - I_80)*100/noI_100

# realative contribution of ICDH with respect to noICDH at 80 % from 80% noICDH 
Diff_rel_80_contrib = M_08_noICDH_abs.subtract(M_08_ICDH_abs, fill_value=0)*100/M_08_noICDH_abs #  (noI_80 - I_80)*100/noI_80



'''
plot data
'''

scale = 0.001
f, axs = plt.subplots(1,3,figsize=(12,3))
plt.subplot(1, 3, 1) 

ax = sns.heatmap(Diff_abs*scale, cmap="Oranges", vmin=0, vmax = 60)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=10)
plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=10)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=10)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_facecolor('lightgray')
ax.tick_params(bottom=False)


plt.subplot(1, 3, 2)
ax =sns.heatmap(Diff_abs_contrib*scale, cmap="Blues",vmin=0, vmax = 2.6) 
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=10)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=10)
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_facecolor('lightgray')
ax.tick_params(left=False,bottom = False, labelleft=False)


plt.subplot(1, 3, 3)
ax = sns.heatmap(Diff_rel_contrib, cmap="Greens",vmin=0, vmax = 0.4)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=10)
plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=10)
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_facecolor('lightgray')
ax.tick_params(left=False,labelleft=False,bottom = False)


plt.savefig('<add path and name of plot file>',dpi=600, bbox_inches='tight')     
