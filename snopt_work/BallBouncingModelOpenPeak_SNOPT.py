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

class BallBouncingModel_Snopt(object):
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
                                    self.num_nodes, endpoint=False, dtype=int)
             
        
    def objective(self, x):
        #
        # The callback for calculating the objective
        # Objective here is the sum of 4th power of external force
        # Objective is scaled down to avoid diverge of optimization
        #
        f_u = sum(x[self.i_u]**4)
        obj = 0.0001*self.interval*f_u
        
        return obj
    
	 
    def gradient(self, x):
        #
        # The callback for calculating the gradient
        # calculate gradiant based on objective function
        #
        grad = np.zeros_like(x)
        grad[self.i_u] = 0.0001*4*self.interval*(x[self.i_u]**3)
        
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
					f, dfdx, dfdxdot, dfdp = self.BallBouncingDynamics(x_a,
															    (x_a - x_p)/h,
                                                                con_a)
            else:
					print 'Do not have the Intergration Method code'
					
            cons[C*p : C*(p+1)] = f

        return cons
    
	
    def jacobianstructureAG(self):
        # Assume the sturcuture of Jac will not change in each iteration
		# (approaved by other project), sparsity functions are used to get
		# the structure of Jac.
		
		# Random initial guess is used to get the jacobian and its structure, 
		# Theoritically, the sturcuture can be get by one time only.
    
      # This code returns jacobian structure of Linear and nonlinear constraints.
      # Since Ipopt requres to separate linear and nonlinear constraints, therefore,
      # Jacobin structures of them need be separate also.
        
        N = self.num_nodes
        h = self.interval
        S = self.num_states
        C = self.num_cons
        U = self.num_u

        Jac_x = np.zeros((C, 2*S))
        
        Jac_c = np.zeros((C, U))
		
        rowA = np.array([])
        colA = np.array([])
        
        rowG = np.array([])
        colG = np.array([])

        for p in range(N-1):
            
            for q in range(5):
                np.random.seed()
                x_p = np.random.random(S)
                np.random.seed()
                x_a = np.random.random(S)
                np.random.seed()
                con_a = np.random.random(U)
    			
                if self.intergration_method == 'backward euler':
                    f, dfdx, dfdxdot, dfdu = self.BallBouncingDynamics(x_a,
                                                                      (x_a - x_p)/h,
                                                                       con_a)
                    
                    Jac_x[:, :S] +=  -dfdxdot/h
                    Jac_x[:, S:2*S] += dfdx + dfdxdot/h
                    Jac_c[:, :U] += dfdu
                    
                else:
                    print 'Do not have the Intergration Method code'
                
            for k in range(C):
                if k < int(C/2):
                    row_xA, col_xA, RA_Jac_xA = find(Jac_x[k, :])
                    row_xAf = row_xA + p*C + k
                    col_xAf = col_xA + p*S
                    
                    row_c, col_c, RA_Jac_c = find(Jac_c[k, :])
                    row_cf = row_c + p*C + k
                    col_cf = col_c + N*S + (p+1)*U
    					
                    rowA = np.hstack((rowA, np.hstack((row_xAf, row_cf))))
                    colA = np.hstack((colA, np.hstack((col_xAf, col_cf))))
                    
                else:
                    row_xG, col_xG, RA_Jac_xG = find(Jac_x[k, :])
                    row_xGf = row_xG + p*C + k
                    col_xGf = col_xG + p*S
                    
                    row_c, col_c, RA_Jac_c = find(Jac_c[k, :])
                    row_cf = row_c + p*C + k
                    col_cf = col_c + N*S + (p+1)*U
    					
                    rowG = np.hstack((rowG, np.hstack((row_xGf, row_cf))))
                    colG = np.hstack((colG, np.hstack((col_xGf, col_cf))))
		
        return rowA, colA, rowG, colG

		
    def jacobianG(self, x):
        #
        # The callback for calculating the nonlinear Jacobian
        # This is nonlinear dynamic functions and first derivatives of Ball Bouncing system
        # 
		
        N = self.num_nodes
        P = self.num_par
        h = self.interval
        S = self.num_states
        C = self.num_cons
        T = self.num_u
        
        con = x[-P:]

        JacG = np.array([])
		
        Jac_p = np.zeros((C, 2*S+T))

        for p in range(N-1):
			
            x_p = x[p*S : (p+1)*S]
            x_a = x[(p+1)*S : (p+2)*S]

            con_a = con[(p+1)*T:(p+2)*T]
			
            if self.intergration_method == 'backward euler':
                f, dfdx, dfdxdot, dfdu = self.BallBouncingDynamics(x_a,
														  (x_a - x_p)/h,
                                                            con_a)
                    
                Jac_p[:, :S] = -dfdxdot/h
                Jac_p[:, S:2*S] = dfdx + dfdxdot/h
                Jac_p[:, 2*S:2*S+T] = dfdu

                for r in range(int(C/2), C):
                    row, col, RA_xG = find(Jac_p[r, :])
                    JacG = np.hstack((JacG, RA_xG))
                        
            else:
				print 'Do not have the Intergration Method code'
				
        return JacG
    
    def jacobianA(self, x):
        #
        # The callback for calculating the linear Jacobian
        # This is linear dynamic functions and first derivatives of Ball Bouncing system
        #
		
        N = self.num_nodes
        P = self.num_par
        h = self.interval
        S = self.num_states
        C = self.num_cons
        T = self.num_u
        
        con = x[-P:]
                
        JacA = np.array([])
		
        Jac_p = np.zeros((C, 2*S+T))

        for p in range(N-1):
			
            x_p = x[p*S : (p+1)*S]
            x_a = x[(p+1)*S : (p+2)*S]

            con_a = con[(p+1)*T:(p+2)*T]
			
            if self.intergration_method == 'backward euler':
                f, dfdx, dfdxdot, dfdu = self.BallBouncingDynamics(x_a,
														  (x_a - x_p)/h,
                                                            con_a)
                    
                Jac_p[:, :S] = -dfdxdot/h
                Jac_p[:, S:2*S] = dfdx + dfdxdot/h
                Jac_p[:, 2*S:2*S+T] = dfdu

                for r in range(0, int(C/2)):
                    row, col, RA_xA = find(Jac_p[r, :])
                    JacA = np.hstack((JacA, RA_xA))
                        
            else:
				print 'Do not have the Intergration Method code'
				
        return JacA
        
    
    def BallBouncingDynamics(self, xs, xsd, u):
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
        
        Fk, dFkdx, dFkdxd = self.SpringDamper(x, xd)
        
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
        

