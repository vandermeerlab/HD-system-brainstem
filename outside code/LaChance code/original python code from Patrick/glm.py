# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:47:16 2024

@author: plachanc
"""

from scipy.sparse import csr_matrix, spdiags
from scipy.stats import wilcoxon
from scipy.optimize import minimize
import numpy as np
import math
import copy
from itertools import chain, combinations


''' here we set various parameter values '''
max_ahv = 180 #deg/s
ahv_bin_size = 18 #deg/s

hd_bin_size = 12 #deg

max_speed = 20 #not sure units on this one
speed_bin_size = 2 #

max_pupil_pos = 30 #deg
pupil_pos_bin_size = 6 #deg

''' calculations based on those '''

n_ahv_bins = int(2*max_ahv/ahv_bin_size)
n_hd_bins = int(360/hd_bin_size)
n_speed_bins = int(max_speed/speed_bin_size)
n_pos_bins = int(2*max_pupil_pos/pupil_pos_bin_size)


def make_X(hds,ahvs,speeds,pupil_pos):
    
    hd_bins = np.digitize(hds,np.linspace(-180,180,int(360/hd_bin_size),endpoint=False)) - 1
    ahv_bins = np.digitize(ahvs,np.linspace(-max_ahv,max_ahv,int(2*max_ahv/ahv_bin_size),endpoint=False)) - 1
    
    speeds[speeds<0] = 0
    speed_bins = np.digitize(speeds,np.linspace(0,max_speed,int(max_speed/speed_bin_size),endpoint=False)) - 1
    
    pupil_pos[pupil_pos<-max_pupil_pos] = -max_pupil_pos
    pupil_pos_bins = np.digitize(pupil_pos,np.linspace(-max_pupil_pos,max_pupil_pos,int(2*max_pupil_pos/pupil_pos_bin_size),endpoint=False)) - 1

    Xhd = np.zeros((len(hds),n_hd_bins))
    Xahv = np.zeros((len(ahvs),n_ahv_bins))
    Xspeed = np.zeros((len(speeds),n_speed_bins))
    Xpos = np.zeros((len(pupil_pos),n_pos_bins))

    for i in range(len(hd_bins)):
        Xhd[i][hd_bins[i]] = 1.
        Xahv[i][ahv_bins[i]] = 1.
        Xspeed[i][speed_bins[i]] = 1.        
        Xpos[i][pupil_pos_bins[i]] = 1.

    X = np.concatenate((Xhd,Xahv,Xspeed,Xpos),axis=1)
    X = csr_matrix(X)

    return X,csr_matrix(Xhd),csr_matrix(Xahv),csr_matrix(Xspeed),csr_matrix(Xpos)

def compute_diags():
    ''' create diagonal matrices for grouped penalization -- implementation 
    modified from Hardcastle 2017 '''
    
    'diagonal matrix for computing differences between adjacent HD bins'
    pos_ones = np.ones(n_hd_bins)
    circ1 = spdiags([-pos_ones,pos_ones],[0,1],n_hd_bins-1,n_hd_bins)
    hd_diag = circ1.T * circ1
    hd_diag=np.asarray(hd_diag.todense())
    hd_diag[0] = np.roll(hd_diag[1],-1)
    hd_diag[n_hd_bins-1] = np.roll(hd_diag[n_hd_bins-2],1)
    
    'one for AHV'
    pos_ones = np.ones(n_ahv_bins)
    ahv1 = spdiags([-pos_ones,pos_ones],[0,1],n_ahv_bins-1,n_ahv_bins)
    ahv_diag = ahv1.T * ahv1
    ahv_diag = np.asarray(ahv_diag.todense())
    
    'and wheel speed'
    pos_ones = np.ones(n_speed_bins)
    speed1 = spdiags([-pos_ones,pos_ones],[0,1],n_speed_bins-1,n_speed_bins)
    speed_diag = speed1.T * speed1
    speed_diag = np.asarray(speed_diag.todense())
    
    'and pupil position'
    pos_ones = np.ones(n_pos_bins)
    pos1 = spdiags([-pos_ones,pos_ones],[0,1],n_pos_bins-1,n_pos_bins)
    pos_diag = pos1.T * pos1
    pos_diag = np.asarray(pos_diag.todense())
        
    return hd_diag, ahv_diag, speed_diag, pos_diag


def objective(params,X,spike_train,smoothers,smooth=True):
    ''' objective function for GLM '''    
    
    u = X * params
    rate = np.exp(u)
    
    f = np.sum(rate - spike_train * u)
    grad = X.T * (rate - spike_train)
    
    if smooth:
        fpen,fgrad = penalize(params,X,spike_train,smoothers)
        f += fpen
        grad += fgrad
    
    # print(f)
    return f,grad

def penalize(params,X,spike_train,smoothers):
    ''' penalize parameter vectors for not being smooth enough '''
    
    hd_diag = smoothers[0]
    ahv_diag = smoothers[1]
    speed_diag = smoothers[2]
    pos_diag = smoothers[3]
    
    hdbeta = 20.
    ahvbeta = 20.
    speedbeta = 20.
    posbeta = 20.
    
    hd_params = params[:n_hd_bins]
    ahv_params = params[n_hd_bins:(n_hd_bins + n_ahv_bins)]
    speed_params = params[(n_hd_bins + n_ahv_bins):(n_hd_bins + n_ahv_bins + n_speed_bins)]
    pos_params = params[(n_hd_bins + n_ahv_bins + n_speed_bins):]
    
    f = np.sum(hdbeta * .5 * np.dot(hd_params.T, hd_diag) * hd_params )
    f += np.sum(ahvbeta * .5 * np.dot(ahv_params.T, ahv_diag) * ahv_params )
    f += np.sum(speedbeta * .5 * np.dot(speed_params.T, speed_diag) * speed_params )
    f += np.sum(posbeta * .5 * np.dot(pos_params.T, pos_diag) * pos_params )

    grad = hdbeta * np.dot(hd_diag,hd_params)
    grad = np.concatenate((grad,ahvbeta * np.dot(ahv_diag,ahv_params)))
    grad = np.concatenate((grad,speedbeta * np.dot(speed_diag,speed_params)))
    grad = np.concatenate((grad,posbeta * np.dot(pos_diag,pos_params)))

    return f,grad


def split_data(X,Xhd,Xahv,Xspeed,Xpos,spike_train,fold):
    ''' split into train and test sets '''

    break_points = np.linspace(0,len(spike_train),51).astype(np.int)


    slices = np.r_[break_points[fold]:break_points[fold + 1],break_points[fold + 10]:break_points[fold + 11],break_points[fold + 20]:break_points[fold + 21],
                          break_points[fold + 30]:break_points[fold + 31],break_points[fold + 40]:break_points[fold + 41]]
    
    test_spikes = spike_train[slices]
    test_Xhd = csr_matrix(Xhd.todense()[slices])
    test_Xahv = csr_matrix(Xahv.todense()[slices])
    test_Xspeed = csr_matrix(Xspeed.todense()[slices])
    test_Xpos = csr_matrix(Xpos.todense()[slices])
    
    train_spikes = np.delete(spike_train,slices,axis=0)
    train_X = csr_matrix(np.delete(X.todense(),slices,axis=0))
    train_Xhd = csr_matrix(np.delete(Xhd.todense(),slices,axis=0))
    train_Xahv = csr_matrix(np.delete(Xahv.todense(),slices,axis=0))
    train_Xspeed = csr_matrix(np.delete(Xspeed.todense(),slices,axis=0))
    train_Xpos = csr_matrix(np.delete(Xpos.todense(),slices,axis=0))
    
    return test_spikes,test_Xhd,test_Xahv,test_Xspeed,test_Xpos,train_spikes,train_X,train_Xhd,train_Xahv,train_Xspeed,train_Xpos
    

def calc_scale_factor(model,hd_params,ahv_params,speed_params,pos_params,train_Xhd,train_Xahv,train_Xspeed,train_Xpos,train_spikes):
    
    u = np.zeros(len(train_spikes))
    
    if 'hd' in model:
        u += train_Xhd * hd_params
    if 'ahv' in model:
        u += train_Xahv * ahv_params
    if 'speed' in model:
        u += train_Xspeed * speed_params
    if 'pos' in model:
        u += train_Xpos * pos_params
    
    rate = np.exp(u)
    
    scale_factor = np.sum(train_spikes)/np.sum(rate)
    
    return scale_factor


def run_final(model,scale_factor,hd_params,ahv_params,speed_params,pos_params,Xhd,Xahv,Xspeed,Xpos,spike_train):
    
    if model != 'uniform':
    
        u = np.zeros(len(spike_train))

        if 'hd' in model:
            u += Xhd * hd_params
        if 'ahv' in model:
            u += Xahv * ahv_params
        if 'speed' in model:
            u += Xspeed * speed_params
        if 'pos' in model:
            u += Xpos * pos_params
        
        rate = np.exp(u) * scale_factor
        
    else:
        
        rate = np.full(len(spike_train),np.mean(spike_train))

    f = -np.sum(rate - spike_train*np.log(rate))
    
    #start array for log-factorials
    lgammas = np.zeros(len(spike_train))
        
    #for every time point...
    for h in range(len(spike_train)):
        #calculate the log-factorial of the number of spikes during that frame,
        #using base 2 log and gamma function on (n_spikes + 1) so numba can work
        lgammas[h] = np.log(math.gamma(spike_train[h]+1))
        
    #subtract the log-factorials
    f -= np.sum(lgammas)
    
    #change from nats to bits
    f = f/np.log(2)

    llps = f/np.sum(spike_train)

    cdict = {}
    #add relevant variables to the appropriate dictionary
    cdict['ll'] = f
    if np.sum(spike_train) > 0:
        cdict['llps'] = float(f/np.sum(spike_train))
    else:
        cdict['llps'] = f
    cdict['lambda'] = rate
    cdict['test_spikes'] = spike_train
    cdict['tot_spikes'] = np.sum(spike_train)
    cdict['scale_factor'] = scale_factor

    return cdict


def get_all_models(variables):
    ''' convenience function for calculating all possible
    combinations of nagivational variables '''
    
    def powerset(variables):
        return list(chain.from_iterable(combinations(variables, r) for r in range(1,len(variables)+1)))
    
    all_models = powerset(variables)
    
    for i in range(len(all_models)):
        all_models[i] = frozenset(all_models[i])
    
    return all_models


def select_model(cdict):
    ''' perform forward search procedure to choose the best model based on
    cross-validation results '''
    
    #variables we'll be looking at
    variables = [('hd',),('ahv',),('speed',),('pos',)]
    #models we'll start with are single-variable models
    models = copy.deepcopy(variables)
    #we haven't found a best model yet so set to NONE
    best_model = None
        
    #while we haven't found a winning model...
    while best_model == None:
        
        #make models frozensets to make things easier
        for i in range(len(models)):
            models[i] = frozenset(models[i])
        
        #start dict for llps measures
        ll_dict = {}
        
        #for each fold in the cross-val...
        for fold in range(10):
            #for each model...
            for modeltype in models:
                #start an entry for that model if we're just starting
                if fold == 0:
                    ll_dict[modeltype] = []
                #collect the llps increase compared to the uniform model
                ll_dict[modeltype].append(cdict[fold][modeltype]['llps']-cdict[fold]['uniform']['llps'])
        
        #make a dict that contains the median llps value for each model
        median_dict = {}
        for modeltype in models:
            median_dict[modeltype] = np.median(ll_dict[modeltype])
        #let's say the potential best new model is the one with the highest median score
        top_model = max(median_dict.keys(), key=(lambda key: median_dict[key]))
        
        # print(top_model)

        #if the top model is a single variable...
        if len(top_model) == 1:
            #set the top model llps data as 'last_model' data
            last_model = ll_dict[top_model]
            #set the top model the 'last_modeltype'
            last_modeltype = top_model
            #create the next set of models -- the current best model plus each new variable
            #-- then start over
            models = []
            for var in variables:
                if var[0] not in top_model:
                    models.append(frozenset(chain(list(top_model),list(var))))
        #otherwise...
        else:
            #use wilcoxon signed ranks to see if the new model is better than the last model
            w,p = wilcoxon(last_model,ll_dict[top_model])
            # print(p)
            #if the new model is better...
            if np.median(last_model) < np.median(median_dict[top_model]) and p < .1:
                #if we can't add any more variables...
                if len(top_model) == len(variables):
                    #test to see if the top model is better than the null model
                    w,p = wilcoxon(ll_dict[top_model],np.zeros(10))
                    # print(p)
                    #if it is, then this is the best model!
                    if np.median(ll_dict[top_model]) > 0 and p < .1:
                        best_model = top_model
                    #otherwise, the cell is unclassifiable
                    else:
                        best_model = 'uniform'
                    
                #otherwise, set this model's llps data as 'last_model' data
                last_model = ll_dict[top_model]
                #set this modeltype as 'last_modeltype'
                last_modeltype = top_model
                #make new set of models -- current best model plus each new variable
                #-- then start over
                models = []
                for var in variables:
                    if var[0] not in top_model:
                        models.append(frozenset(chain(list(top_model),list(var))))
            #otherwise, the best model is probably the last model
            else:
                #check if the last model is better than the null  model
                w,p = wilcoxon(last_model,np.zeros(10))
                # print(p)
                # print(last_model)
                #if it is, then this is the best model!
                if np.median(last_model) > 0 and p < .1:
                    best_model = last_modeltype
                #otherwise, the cell is unclassifiable
                else:
                    best_model = 'uniform'
                    
    #return the best model
    return best_model


def run_classifier(hds,ahvs,speeds,pupil_pos,spike_train):
    ''' run 10-fold cross-validation '''
    
    variables = [('hd'),('ahv'),('speed'),('pos')]
    all_models = get_all_models(variables)
    
    #make design matrix for basic model
    X,Xhd,Xahv,Xspeed,Xpos = make_X(hds,ahvs,speeds,pupil_pos)

    #make smoothing matrices for regularization in cross-validation
    smoothers = compute_diags()
                
    cdict = {}

    for fold in range(10):
        cdict[fold] = {}
        test_spikes,test_Xhd,test_Xahv,test_Xspeed,test_Xpos,train_spikes,train_X,train_Xhd,train_Xahv,train_Xspeed,train_Xpos = split_data(X,Xhd,Xahv,Xspeed,Xpos,spike_train,fold)
        params = np.zeros(np.shape(train_X)[1])
        result = minimize(objective,params,args=(train_X,train_spikes,smoothers),jac=True,method='L-BFGS-B')

        params = result.x
        
        hd_params = params[:n_hd_bins]
        ahv_params = params[n_hd_bins:(n_hd_bins + n_ahv_bins)]
        speed_params = params[(n_hd_bins + n_ahv_bins):(n_hd_bins + n_ahv_bins + n_speed_bins)]
        pos_params = params[(n_hd_bins + n_ahv_bins + n_speed_bins):]
        
        for model in all_models:
            
            scale_factor = calc_scale_factor(model,hd_params,ahv_params,speed_params,pos_params,train_Xhd,train_Xahv,train_Xspeed,train_Xpos,train_spikes)
            cdict[fold][model] = run_final(model,scale_factor,hd_params,ahv_params,speed_params,pos_params,test_Xhd,test_Xahv,test_Xspeed,test_Xpos,test_spikes)
        
        cdict[fold]['uniform'] = run_final('uniform',1.,hd_params,ahv_params,speed_params,pos_params,test_Xhd,test_Xahv,test_Xspeed,test_Xpos,test_spikes)

    best_model = select_model(cdict)
    
    return best_model
