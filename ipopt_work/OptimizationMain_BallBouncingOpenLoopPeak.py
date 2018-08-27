#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 17:59:21 2017

Revised on Wed May 02 12:13:00 2018

This is the main code to optimize Ball Bouncing problem.

Change num_nodes and duration can optimize different length of bouncing trajectory.

Usually, inverval is keeped as 10 ms, therefore, num_nodes = 100*duration + 1

@author: huawei
"""

import ipopt
import numpy as np
from BallBouncingModelOpenPeak import BallBouncingModel
import time
import os

# define optimization parapmeters 
num_nodes = 1001  # number of direct collocation nodes
duration = 10.0  # duration time of optimization trajectory
dint = int(duration)
interval = duration/(num_nodes - 1)  # calcuate interval between direct collocation
                                     # nodes

num_par = num_nodes # number of control nodes, equal number of nodes because of open control
num_states = 2  # number of states of each node
num_cons = 2  # number of constraints of each node

Scaling = 100.0  # scaling of constraints

# define peak targets and bottom targets
num_peak = int(duration/2.0) + 1  # total number of peak 
peak_dis_ind = np.linspace(0, (num_nodes-1)*2, num_peak, dtype=np.int32) # index of peak targets # 0, 2, 4...
peak_dis = np.linspace(10.0, 10.0, num_peak)  # initialize peak targets with 10
np.random.seed(seed=88)   # define random seeds
peak_dis[1:] = peak_dis[1:] - 2*np.random.random(num_peak-1)  # random chose peak targets, except the first
                                                             # first will keep with 10
bot_dis_ind = np.linspace(200, (num_nodes-1)*2-200, num_peak-1, dtype=np.int32)  # index of bottom peak

# solver options, choose one
#Solver = 'MA27'
#Solver = 'MA57'
#Solver = 'MA77'
#Solver = 'MA86'
Solver = 'MUMPS'
#Solver = 'Pardiso'

i_x = np.linspace(0, num_states*num_nodes, num_nodes,  # x index
                  endpoint=False, dtype=np.int32) + 0  
i_v = np.linspace(0, num_states*num_nodes, num_nodes,  # xdot index
                  endpoint=False, dtype=np.int32) + 1
                                     
i_u = np.linspace(num_states*num_nodes, num_states*num_nodes+num_par,  # controller index
                  num_nodes, endpoint=False, dtype=np.int32) + 0

lb_x = np.zeros((num_states*num_nodes))  # initialize low bounds of optimization
                                         # parameters, include x, xdot, and u

lb_x[i_x] = -0.1   # set lower bounds of x to -0.1 m
lb_x[i_v] = -20.0  # set lower bounds of xdot to -20.0 m/s
lb_x[peak_dis_ind] = peak_dis  # set lower bounds of x at peak target index to previous generated peak numbers
lb_x[peak_dis_ind+1] = 0.0   # set lower bounds of xdot at peak target index to 0 m/s
lb_x[bot_dis_ind] = -0.02  # set lower bounds of x at bottom target index to -0.02 m
lb_c = np.zeros((num_par)) -200.0  # set lower bounds of controller to -100 N
                  
ub_x = np.zeros((num_states*num_nodes))  # initialize upper bounds of optimization
                                         # parameters, include x, xdot, and u
ub_x[i_x] = 10.0  # set upper bounds of x to 10 m
ub_x[i_v] = 20.0  # set upper bounds of xdot to 20.0 m/s
ub_x[peak_dis_ind] = peak_dis  # set upper bounds of x at peak target index to previous generated peak numbers
ub_x[peak_dis_ind+1] = 0.0   # set upper bounds of xdot at peak target index to 0 m/s

ub_x[bot_dis_ind] = -0.02  # set upper bounds of x at bottom target index to -0.02 m
               
ub_c = np.zeros((num_par)) + 200.0  # set upper bounds of controller to 100 N
    
lb = list(np.hstack((lb_x, lb_c)))  # combine states and controller low bounds
ub = list(np.hstack((ub_x, ub_c)))  # combine states and controller upper bounds

cl = list(np.zeros((num_cons*(num_nodes -1), 1)))  # set lower  and bounds of constraints to 0
cu = list(np.zeros((num_cons*(num_nodes -1), 1)))  # since they are all equality constraints

# define ipopt problem according to cyipopt interface settings

nlp = ipopt.problem(
            n=num_states*num_nodes + num_par,  # total number of optimization parameters
            m=num_cons*(num_nodes -1),  # number of total constraints
            
            # optimization problems, include objectives, objective gradiants, constraints, and jacobins
            problem_obj=BallBouncingModel(num_nodes, num_states, num_par, interval,
                                           scaling =Scaling, integration_method='backward euler'),
                                           
            lb=lb,  # lower bounds of optimization parameters
            ub=ub,  # upper bounds of optimization parameters
            cl=cl,  # lower bounds of optimization constraints
            cu=cu   # lower bounds of optimization constraints
            )

# set options of ipopt optimization problems
nlp.addOption('linear_solver', Solver)  # define solver
nlp.addOption('max_iter', 5000)  # define maximum iterations
nlp.addOption('hessian_approximation', 'limited-memory')  # define hessian approximation
                                                         # as limited-memory, since no hassian
                                                        # will be provided
nlp.addOption('tol', 1e-4)  # define tolerance of constraints 
nlp.addOption('acceptable_tol', 1e-3)  # define acceptable tolerance of constraints 
nlp.addOption('max_cpu_time', 4e+4)   # define maximum cup time to solve the problem

#  initialize state parameters with x = 5.0, xdot = 0, x[peak] = oeak
x_init = np.zeros(num_states*num_nodes)
x_init[i_x] = 5.0
x_init[i_v] = 0.0 
x_init[peak_dis_ind] = peak_dis
x_init[bot_dis_ind] = -0.02

unknown = np.zeros(num_par) + 10.3  # initialize controller with 10.3N
x0 = np.hstack((x_init, unknown))  # combine states and controller initialization

Start_Time = time.time()  # record current
x, info = nlp.solve(x0)  # run the optimization with initial guess
Opt_Time = time.time() - Start_Time  # calcuate time consuming of optimization 

Sta = info['status']  # get status of optimization when finish
R_cst = info['g']  # get residue of constraints when finish

# save optimizated optimization parameters, residue of constraints, optimized controller,
# and status/optimization_time into text files.

newpath = 'ResultBallDropPeak' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

FitName = newpath+'/TrajectoryResult_'+Solver+'_'+str(dint)+'_Sca100.txt'
with open(FitName,'w') as Outfile:
    for m in range(0,num_nodes):
        StringP = ""
        for n in range(0,num_states):
            StringP += str(x[m*num_states + n])
            StringP += " "
        StringP += "\n"
        Outfile.write(StringP)
        
RcstName = newpath+'/ResiduleConstraints_'+Solver+'_'+str(dint)+'_Sca100.txt'
with open(RcstName,'w') as Outfile:
    StringP = ""
    for m in range(0,len(cl)):
        StringP += str(R_cst[m])
        StringP += "\n"
    Outfile.write(StringP)

ContrName = newpath+'/Torques_'+Solver+'_'+str(dint)+'_Sca100.txt'
with open(ContrName,'w') as Outfile:
    StringP = ""
    for i in range(0,num_nodes):
        StringP += str(x[num_nodes*num_states + i])
        StringP += "\n"
    Outfile.write(StringP)
    
ContrName = newpath+'/RMS_'+Solver+'_'+str(dint)+'_Sca100.txt'
with open(ContrName,'w') as Outfile:
    StringP = ""
    StringP += str(Sta)
    StringP += "\n"
    StringP += str(Opt_Time)
    Outfile.write(StringP)