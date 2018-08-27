#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 10:46:43 2017

Revised on Wed May 02 12:13:00 2018

@author: huawei

This code defined Ball Bouncing dynamics for optimization with direct collocation.

This code will provide Objective, Gradient, Constraint, Jacobin, and Jacobin strucuture
"""

import numpy as np
from scipy.sparse import find

class BallBouncingModel(object):
    def __init__(self, num_nodes, num_states, num_par, interval, scaling =1, integration_method='backward euler' ):
        
        self.num_nodes = num_nodes
        self.num_states = num_states
        self.num_cons = num_states
        self.num_par = num_par
        self.num_conspernode = num_states
        
        self.num_u = 1
		
        self.interval = interval
        self.scaling = scaling
        self.intergration_method = integration_method
        
        # define index of external force in optimization parameters
        self.i_u = np.linspace(self.num_states*self.num_nodes,
                                   self.num_states*self.num_nodes + self.num_par,
                                    self.num_nodes, endpoint=False, dtype=np.int32)
             
        
    def objective(self, x):
        #
        # The callback for calculating the objective
        # Objective here is the sum of 4th power of external force
        # Objective is scaled down to avoid diverge of optimization
        #
        
        f_u = sum(x[self.i_u]**4)
        obj = 0.0001*self.interval*f_u/self.num_nodes
        
        return obj
    
	 
    def gradient(self, x):
        #
        # The callback for calculating the gradient
        # calculate gradiant based on objective function
        #
        
        grad = np.zeros_like(x)
        grad[self.i_u] = 0.0001*4*self.interval*(x[self.i_u]**3)/self.num_nodes
        
        return grad

		
    def constraints(self, x):
        #
        # The callback for calculating the constraints
        # Calcaute residule of constraint functions based on optimizing parameters
        # In direct collocation here, constraints are system dyanmics functions 
        #
        N = self.num_nodes
        S = self.num_states
        C = self.num_cons
        h = self.interval
        T = self.num_u

        cons = np.zeros(C*(N - 1))
        con = x[-self.num_par:]
        
        for p in range(N-1):
			
            x_p = x[p*S : (p+1)*S]
            x_a = x[(p+1)*S : (p+2)*S]

            con_a = con[(p+1)*T : (p+2)*T]
				
            if self.intergration_method == 'backward euler':
					f, dfdx, dfdxdot, dfdp = self.BallBouncing(x_a,
            													 (x_a - x_p)/h,
                                                                con_a)
            else:
					print 'Do not have the Intergration Method code'
					
            cons[C*p : C*(p+1)] = f

        return cons
    
	
    def jacobianstructure(self):
        # Assume the sturcuture of Jac will not change in each iteration
		# (approaved by other project), sparsity functions are used to get
		# the structure of Jac.
		
		# random initial guess is used to get the jacobian and its structure, 
		# Theoritically, the sturcuture can be get by one time only.
		
        N = self.num_nodes
        P = self.num_par
        h = self.interval
        S = self.num_states
        C = self.num_cons
        U = self.num_u
		
        np.random.seed()
        x = np.random.random(N*S+P)
		
        con = x[-P:]
		
        Jac_x = np.zeros((C, 2*S))
        
        Jac_c = np.zeros((C, U))
		
        row = np.array([])
        col = np.array([])

        for p in range(N-1):
            x_p = x[p*S : (p+1)*S]
            x_a = x[(p+1)*S : (p+2)*S]

            con_a = con[(p+1)*U : (p+2)*U]
			
            if self.intergration_method == 'backward euler':
                f, dfdx, dfdxdot, dfdu = self.BallBouncing(x_a,
														 (x_a - x_p)/h,
                                                            con_a)
                
                Jac_x[:, :S] =  -dfdxdot/h
                Jac_x[:, S:2*S] = dfdx + dfdxdot/h
                Jac_c[:, :U] = dfdu
                
                for k in range(C):
                    row_x, col_x, RA_Jac_x = find(Jac_x[k, :])
                    row_c, col_c, RA_Jac_c = find(Jac_c[k, :])
                    
                    row_xf = row_x + p*C + k
                    row_cf = row_c + p*C + k
                                    
                    col_xf = col_x + p*S
                    col_cf = col_c + N*S + (p+1)*U
					
                    row = np.hstack((row, np.hstack((row_xf, row_cf))))
                    col = np.hstack((col, np.hstack((col_xf, col_cf))))
            else:
				print 'Do not have the Intergration Method code'
                
        row = row.astype(int)
        col = col.astype(int)
        
        return (row, col)

		
    def jacobian(self, x):
        #
        # The callback for calculating the Jacobian
        #
		
        N = self.num_nodes
        P = self.num_par
        h = self.interval
        S = self.num_states
        C = self.num_cons
        T = self.num_u
        
        con = x[-P:]
                
        Jac = np.array([])
		
        Jac_p = np.zeros((C, 2*S+T))

        for p in range(N-1):
			
            x_p = x[p*S : (p+1)*S]
            x_a = x[(p+1)*S : (p+2)*S]

            con_a = con[(p+1)*T:(p+2)*T]
			
            if self.intergration_method == 'backward euler':
                f, dfdx, dfdxdot, dfdu = self.BallBouncing(x_a,
														  (x_a - x_p)/h,
                                                             con_a)
                    
                Jac_p[:, :S] = -dfdxdot/h
                Jac_p[:, S:2*S] = dfdx + dfdxdot/h
                Jac_p[:, 2*S:2*S+T] = dfdu

                for r in range(C):
                    row, col, RA_x = find(Jac_p[r, :])
                    Jac = np.hstack((Jac, RA_x))
            else:
				print 'Do not have the Intergration Method code'
				
        return Jac
        
    
    def BallBouncing(self, xs, xsd, u):
        #
        # This is dynamic functions and first derivatives of Ball Bouncing system
        #
        
        M = 1.0        
        g = 9.8
        
        scaling = self.scaling
        S = self.num_states
        
        x = xs[:S/2]
        xd = xs[S/2:]
        xdv = xsd[:S/2]
        xdd = xsd[S/2:]
        
        Fk, dFkdx, dFkdxd = self.ContactModel(x, xd)
        
        f = np.zeros(S)
        f[0] = xd - xdv
        f[1] = (-M*g - M*xdd + Fk + u)/self.scaling
         
        dfdx = np.zeros((S, S))
        dfdx[0, 0] = 0
        dfdx[0, 1] = 1.0
        dfdx[1, 0] = dFkdx/self.scaling
        dfdx[1, 1] = dFkdxd/self.scaling
            
        dfdxd = np.zeros((S, S))
        dfdxd[0, 0] = -1.0
        dfdxd[0, 1] = 0
        dfdxd[1, 0] = 0
        dfdxd[1, 1] = -M/self.scaling
             
        dfdu = np.zeros((S, 1))
        dfdu[1, 0] = 1.0/scaling
        
        return f, dfdx, dfdxd, dfdu
    
    def ContactModel(self, x, xd):
        #
        # This is the strong nonlinear contact model, no damper was included
        #
        Ka = 0.1
        Kg = 5e7
        C = 0.0
        
        d = 0.5*(abs(x) - x)
        Fk = (Kg*d**3 - Ka*x) - C*xd
        
        if x >=0:
            dFkdx = -Ka
        else:
            dFkdx = (-3*Kg*d**2 - Ka)
            
        dFkdxd = -C
        
        return Fk, dFkdx, dFkdxd
        
    def intermediate(
        self,
        alg_mod,
        iter_count,
        obj_value,
        inf_pr,
        inf_du,
        mu,
        d_norm,
        regularization_size,
        alpha_du,
        alpha_pr,
        ls_trials
        ):
        
        #
        # Example for the use of the intermediate callback.
        #
        print "Objective value at iteration #%d is - %g" % (iter_count, obj_value)
