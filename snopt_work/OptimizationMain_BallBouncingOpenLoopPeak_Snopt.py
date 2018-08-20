"""
This is the main code to optimize Ball Bouncing problem.

Change num_nodes and duration can optimize different length of bouncing trajectory.

Usually, inverval is keeped as 10 ms, therefore, num_nodes = 100*duration + 1

"""

import numpy        as np
from   optimize import snopta, SNOPT_options
from   BallBouncingModelOpenPeak_SNOPT import BallBouncingModel_Snopt
import time
import os

# define optimization parapmeters 
num_nodes = 1001  # number of direct collocation nodes
duration = 10.0  # duration time of optimization trajectory
dint = int(duration)
interval = duration/(num_nodes - 1)  # calcuate interval between direct collocation
                                     # nodes

num_par = num_nodes  # number of control nodes, equal number of nodes because of open control
num_states = 2  # number of states of each node
num_cons = 2  # number of constraints of each node

Scaling = 100.0  # scaling of constraints

# define peak targets and bottom targets
num_peak = int(duration/2.0) + 1  # total number of peak 
peak_dis_ind = np.linspace(0, (num_nodes-1)*2, num_peak, dtype=np.int32)  # index of peak targets # 0, 2, 4...
peak_dis = np.linspace(10.0, 10.0, num_peak)  # initialize peak targets with 10
np.random.seed(seed=88)  # define random seeds
peak_dis[1:] = peak_dis[1:] - 2*np.random.random(num_peak-1)  # random chose peak targets, except the first
                                                              # first will keep with 10

bot_dis_ind = np.linspace(200, (num_nodes-1)*2-200, num_peak-1, dtype=np.int32)  # index of bottom peak

# load BallBouncing model to get objective, gradient, constraints and jacobians
problem=BallBouncingModel_Snopt(num_nodes, num_states, num_par, interval,
                               scaling =Scaling, integration_method='backward euler')


i_x = np.linspace(0, num_states*num_nodes, num_nodes, endpoint=False, dtype=np.int32) + 0  # x index
i_v = np.linspace(0, num_states*num_nodes, num_nodes, endpoint=False, dtype=np.int32) + 1  # xdot index
                                     
i_u = np.linspace(num_states*num_nodes, num_states*num_nodes+num_par, num_nodes, endpoint=False, dtype=np.int32) + 0  # controller index

def BallBouncing_objFG(status,x,needF,F,needG,G):
    #
    # this function combine obejctive and nonlinear constraints, and put objective
    # at the first row of nonlinear constraints.
    # gradiants are combined with nonlinear jacobians also.
    
    F[0] = problem.objective(x)         # objective row
    F[1:(num_nodes-1)*num_cons+1] = problem.constraints(x)
    
    G = np.array([])
    grad = problem.gradient(x)
    G = np.hstack((G, grad[2*num_nodes:]))
    Gg = problem.jacobianG(x)
    G = np.hstack((G, Gg))

    return status, F, G

# index of objective in nonlinear jacobins.
row_obj = np.linspace(0, 0, num_nodes, dtype=np.int32)
col_obj = np.linspace(0, num_nodes-1, num_nodes, dtype=np.int32) + 2*num_nodes

# linear and nonlinear jacobin structures
Arow, Acol, rowG, colG = problem.jacobianstructureAG()

# Full nonlinear jacobin structure 
Grow = np.hstack((row_obj, rowG+1))
Gcol = np.hstack((col_obj, colG))
Gind = ((Grow), (Gcol))

# define a huge number to infinite
inf   = 1.0e10

# set optimization options
options = SNOPT_options()
options.setOption('Verbose',True)
options.setOption('Solution print',True)
options.setOption('Print filename','sntoya_'+str(dint)+'IpoptInit.out')
options.setOption('Summary frequency',1)

# set initial guess
x0      = np.zeros(3*num_nodes)
x0[i_x] = 5.0
x0[i_v] = 0.0
x0[peak_dis_ind] = peak_dis
x0[bot_dis_ind] = -0.02

x0[i_u] = 10.3

# set lower bounds of optimization parameters
xlow    = np.zeros(3*num_nodes)
xlow[i_x] = -0.1
xlow[i_v] = -20.0
xlow[i_u] = -100.0

xlow[peak_dis_ind] = peak_dis
xlow[peak_dis_ind+1] = 0
xlow[bot_dis_ind] = -0.02

# set upper bounds of optimization parameters
xupp    = np.zeros(3*num_nodes)
xupp[i_x] = 10.0
xupp[i_v] = 20.0
xupp[i_u] = 100.0

xupp[peak_dis_ind] = peak_dis
xupp[peak_dis_ind+1] = 0
xupp[bot_dis_ind] = -0.02

# set bounds of constraint functions, zero except the first node (objective)
Flow    = np.zeros((num_nodes-1)*num_cons+1)
Flow[0] = 0
Fupp    = np.zeros((num_nodes-1)*num_cons+1)
Fupp[0] = inf

# optimization information, number of optimization parameters, number of constraints
n       = 3*num_nodes
nF      = (num_nodes-1)*num_cons+1


ObjRow  = 1  # objective function is the first row
Ag = problem.jacobianA(x0)  # Jacobin of linear constraints

A = ((Ag), (Arow+1), (Acol))  # Linear jacobin and structure

# do optimization and count time
rescent_time = time.time()
result = snopta(BallBouncing_objFG,n,nF,x0=x0,name='sntoyaFG_'+str(dint)+'IpoptInit',xlow=xlow,xupp=xupp,
                Flow=Flow,Fupp=Fupp,ObjRow=ObjRow, A = A, G = Gind)
Time_sol = time.time() - rescent_time

newpath = 'ResultBallDropPeak' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

FitName = newpath+'/TrajectoryResult_'+str(dint)+'_Sca100_Preinit.txt'
with open(FitName,'w') as Outfile:
    for m in range(0,num_nodes):
        StringP = ""
        for n in range(0,num_states):
            StringP += str(result.x[m*num_states + n])
            StringP += " "
        StringP += "\n"
        Outfile.write(StringP)
        
RcstName = newpath+'/ResiduleConstraints_'+str(dint)+'_Sca100_Preinit.txt'
with open(RcstName,'w') as Outfile:
    StringP = ""
    for m in range(0,len(result.F)-1):
        StringP += str(result.F[m+1])
        StringP += "\n"
    Outfile.write(StringP)

ContrName = newpath+'/Torques_'+str(dint)+'_Sca100_Preinit.txt'
with open(ContrName,'w') as Outfile:
    StringP = ""
    for i in range(0,num_nodes):
        StringP += str(result.x[num_nodes*num_states + i])
        StringP += "\n"
    Outfile.write(StringP)
    
ContrName = newpath+'/RMS_'+str(dint)+'_Sca100_Preinit.txt'
with open(ContrName,'w') as Outfile:
    StringP = ""
    StringP += str(result.info)
    StringP += "\n"
    StringP += str(result.objective)
    StringP += "\n"
    StringP += str(Time_sol)
    StringP += "\n"
    StringP += str(result.major_itns)
    StringP += "\n"
    StringP += str(result.iterations)
    Outfile.write(StringP)