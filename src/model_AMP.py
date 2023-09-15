'''

name:   model_abiotic_batch.py 

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src'

author: DKM

goal: Loop model of Monoculture BCC assays to graph 0 H phyotplankton biomass and model of said biomass via odelib

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
plt.rcParams.update({'font.size': 16})
plt.rcParams['legend.fontsize'] = 'small'


#####################################################
# read in data and formatting for modeling and graphing
#####################################################

#main df read in 
df_all = pd.read_csv("../data/AMPMedia.csv",header=0)


#format empty columns and column names 
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time (hr)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#df_all[df_all['log1','log2','log3']] = np.log(df_all[df_all['rep1','rep2','rep3']])

#making log of data to look at error
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)

exps = df_all['ID'].unique()
nexps = exps.shape[0]   #making numnber list for exp



#####################################################
#functions  for modeling and graphing model uncertainty 
#####################################################
def set_best_params(model,posteriors,snames):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    bestchain = posteriors.iloc[im]["chain#"]
    posteriors = posteriors[posteriors["chain#"]==bestchain]
    model.set_parameters(**posteriors.loc[im][a0.get_pnames()].to_dict())
    model.set_inits(**{o:posteriors.loc[im][a0.get_pnames()].to_dict()[o+'0'] for o in ['H']})
###############
#####only set for 0 a for idk if 400 model is working correctly. #######

#function for plotting uncertainty once model has been run 
def plot_uncertainty(ax,model,posteriors,ntimes):
        for a in range(ntimes):
            im = rd.choice(posteriors.index) 
            model.set_inits(**{'H':posteriors.loc[im][model.get_pnames()].to_dict()['H0']})
            model.set_parameters(**posteriors.loc[im][model.get_pnames()].to_dict())
            mod = model.integrate()
            ax.plot(mod.time,mod['H'],c=str(0.8),lw=1,zorder=1)


#actual model that will be run by the odelib model framework
def abiotic(y,t,params):
    deltah,Sh = params[0], params[1]
    H = y[0]
    dHdt = Sh - deltah*H 
    return [dHdt]



#initiating the model as a class in odelib (give us use of the methods in this class - like integrate :-) 
def get_model(df):
    a1=ODElib.ModelFramework(ODE=abiotic,
                             parameter_names=['deltah','Sh', 'H0'],
                             state_names = snames,
                             dataframe=df,
                             deltah = deltah_prior.copy(),
                             Sh = Sh_prior.copy(),
                             H0  = H0_prior.copy(),
                             t_steps=1000,
                             H = H0_mean,
                             )
    return a1

 
#find closesst time 
def get_residuals(self):
    mod = self.integrate(predict_obs=True)
    res = (mod.abundance - self.df.abundance)   #this is not same species 
    mod['res'] = res
    return(mod)


#set up figures to populate from loop 

fig3,ax3 = plt.subplots(nexps,3,figsize=[20,14]) #plot creation and config 
fig3.suptitle('Abiotic HOOH Model for AMP-A') #full title config


#####################################################
# Set up large loop for going through all exps in DF
#####################################################

for (e,ne) in zip(exps,range(nexps)):  #looping through each exp with number 
    df0 = df_all[(df_all['ID'] == e)] #setting working df as a single Experiment in df_all
    inits0 = pd.read_csv("../data/inits/AMP0.csv") #use and updat inits 

#for each exp this is done so the inits can update within the fitting loop

#####################################################
#model param and state variable set up 
#####################################################
    nits = 10000 #number of iterations for MCMC (increase number for more bell curve of param distributions)

    snames = ['H']  # state variable names
    pw = 1 #sigma we give model parameter]

    Hstart = df0.loc[df0['time'] == 0, 'abundance'].iloc[0] #Setting Hscale start as initial condition of each DF
    deltah_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,\
                                  hyperparameters={'s':pw,'scale':0.2}) #setting prior gruess and sigma of deltah param
    Sh_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm, 
                              hyperparameters={'s':pw,'scale':1}) #setting prior gruess and sigma of Sh param
    H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,\
                              hyperparameters={'s':(pw/5),'scale':Hstart}) #setting state variiable  prior guess
    H0_mean = 2 #setting H mean for odelib search 

#####################################################
# initialize and run model 
#####################################################

    a0 = get_model(df0) #initialize model from df 

# do fitting (get posteriors) for model 
    posteriors0 = a0.MCMC(chain_inits=inits0,iterations_per_chain=nits,cpu_cores=1, print_report=False) #generating optimized posteriors for multiple param choices
    set_best_params(a0,posteriors0,snames) # set best params for model using function difined above  #getting params that give lowest error between model posterior set and data
    mod0 = a0.integrate() # run model with optimal params for 0 and save in mod

#get residuals from model 
    a0res = get_residuals(a0)  #is this using the best fit or just a first run???

#########################################################
# graphing df and model outputs together
#########################################################

# Set up graph for Dynamics and param histograms
    fig1,ax1 = plt.subplots(2,1,sharex=True, figsize=[8,5])  #making fig for param v param graphing 
    fig1.suptitle('deltah vs Sh in Exp'+ str(e))
    ax1[0].set_ylabel('HOOH Concentration nM/mL')
    ax1[0].set_xlabel(' Time (hrs) ')
    ax1[0].semilogy()
    l1 = ax1[0].legend(loc = 'best')
    l1.draw_frame(False)
    ax1[1].set_xlabel('Frequency Sh')
    ax1[0].set_ylabel('H concentration')
    ax1[1].set_ylabel('Frequency deltah')

    #graphing each assay's parameters against each other 
    ax1[0].plot(df0.time,df0.abundance, marker='o',label = 'AMP data ' + str(e) )
    ax1[0].plot(mod0.time,mod0['H'],c='r',lw=1.5,label=' model best fit')
    ax1[1].scatter((posteriors0.Sh),(posteriors0.deltah))
    fig1.savefig('../figures/AMP_'+str(e)+'_params.png')
    
    #################################
    #graphing logged parameter values
    ##################################
    #crating and config of fig 
    fig2,ax2 = plt.subplots(2,1,sharex=True,figsize=[8,5]) #make plot
    fig2.suptitle('Trace plots for '+ str(e)) #set main title 
    fig2.subplots_adjust(right=0.90, wspace = 0.25, top = 0.85) #shift white space for better fig view
    fig2.supxlabel('Model Iteration') #set overall x title 
    ax2[0].set_title('deltah')
    ax2[0].set_ylabel('log deltah')
    ax2[1].set_title('Sh ')
    ax2[1].set_ylabel('log Sh')
    #graphing iteration number vs parameter numbert logged 
    ax2[1].scatter(posteriors0.iteration,np.log(posteriors0.Sh))
    ax2[0].scatter(posteriors0.iteration,np.log(posteriors0.deltah))
    
    fig2.savefig('../figures/AMP_'+str(e)+'_Trace.png')
    
    #########################################
    #graphing Residuals of best model vs data 
    ##########################################

    #making and confing of residuals plot
    fig4,ax4 = plt.subplots(2,1,sharex = True,figsize=[8,5])
    fig4.suptitle('Residuals vs Fit Value ')
    fig4.supylabel('Model Value (H)')
    fig4.supxlabel('Residual')

    #config legends for data differentialtion 
    l4 = ax4[0].legend()
    l5 = ax4[1].legend()
    l4.draw_frame(False)
    l5.draw_frame(False)

    #plotting residual function output residual and abundance columns 
    ax4[0].scatter(a0res['res'], a0res['abundance'],label = '0 H') #where )
    fig4.savefig('../figures/AMP_'+str(e)+'_residuals.png')

    #graph large fig of all runs together 
    ax3[ne,0].plot(df0.time,df0.abundance, marker='o',label = 'AMP data ' + str(e) ) #data of 0 H assay
    ax3[ne,0].plot(mod0.time,mod0['H'],c='r',lw=1.5,label=' Model best fit') #best model fit of 0 H assay
    plot_uncertainty(ax3[ne,0],a0,posteriors0,100) #plotting 100 itterations of model search for 0 H assay 
    ax3[2,0].set_ylabel('HOOH Concentration nM/mL')
    ax3[(ne-1),0].set_xlabel(' Time (hrs) ')
    ax3[ne,0].semilogy()
    l3 = ax3[ne,0].legend(loc = 'upper left')
    l3.draw_frame(False)
    ax3[ne,1].hist((np.log(posteriors0.Sh)))
    ax3[ne,2].hist((np.log(posteriors0.deltah)))
    ax3[0,1].set_title('Sh')
    ax3[0,2].set_title('deltah')
    
plt.show()


fig4.savefig('../figures/AMP_all_dynamics.png')



# 'program finished' flag
print('\n Done my guy \n')

print('\n Im free Im free! Im done calculating!' )



