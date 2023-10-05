'''

name:   vis_MHMpy 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: abiotic media producton of heterotrophic media ()

working on: ln of data in df for uncertainty, loop of all dfs in df_all for model and intits? 

'''

#read in needed packages 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ODElib
import random as rd


#####################################################
#set figure RC params 
#####################################################
plt.rcParams["figure.dpi"] = 300
plt.rcParams.update({'font.size': 18})
plt.rcParams['legend.fontsize'] = 'small'


#####################################################
# read in data and formatting for modeling and graphing
#####################################################

#main df read in 
df_all = pd.read_csv("../data/MHM_blank.csv",header=0)


#format empty columns and column names 
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time (hr)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#df_all[df_all['log1','log2','log3']] = np.log(df_all[df_all['rep1','rep2','rep3']])

#making log of data to look at error
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['sigma'] = np.std(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)

exps = df_all['ID'].unique()
exps = exps[~pd.isna(exps)]

nexps = exps.shape[0]   #making numnber list for exp


#####################################################
#set up figures to populate from loop 
#####################################################
#####################################################
fig1,ax1 = plt.subplots(nexps,3,figsize=[15,12]) #plot creation and config 
fig1.suptitle('MHM Data Dynamics') #full title config
ax1[0,0].set_title('HOOH Production ')


fig1.subplots_adjust(left=0.15, bottom=0.10, right=0.90, top=0.9, wspace=0.30, hspace=0.30) #shift white space for better fig view


#####################################################
# Set up large loop 
#####################################################

for (e,ne) in zip(exps,range(nexps)):  #looping through each exp with number 
    df = df_all[(df_all['ID'] == e)] #setting working df as a single Experiment in df_all
    ax1[ne,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker = '*', c='purple',label =  str(e) + ' mean' ) #data of 0 H assay
    ax1[ne,0].errorbar(df.time,df.log_abundance, yerr=df.log_sigma, marker = '*', c='pink',label =  str(e) + ' log mean' )
    ax1[ne,1].scatter(df.abundance,df.sigma, c='k')
    ax1[ne,2].scatter(df.log_abundance,df.log_sigma, c='b')
    ax1[ne,0].semilogy()
    l1 = ax1[ne,0].legend(loc = 'center right')
    l1.draw_frame(False)
    

for a in ax1[:,0]:
    a.set_ylabel('HOOH nM')
    a.set_xlabel('Time (days)')

# ylabels
for a in ax1[:,1]:
    a.set_ylabel('Abundance')
    a.set_xlabel('STDV')

for a in ax1[:,2]:
    a.set_ylabel('Abundance')
    a.set_xlabel('STDV')



plt.show()


fig1.savefig('../figures/MHM_data_dynamics.png')




# 'program finished' flag
print('\n ~~~****~~~****~~~ \n')
print(' Done my guy ')
print('\n ~~~****~~~****~~~ \n')
print('\n Im free Im free! Im done calculating!' )
print('\n ~~~****~~~****~~~ \n')


