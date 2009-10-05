from math import *
import numpy as np
import scipy as sc
import scipy.linalg
import pylab as pylab
from scipy import *
from pylab import *
from Sjoe_m_func2 import *
from Sjoe_index import *
from qrp import qrp
import IVP
from numpy import max, rank

np.set_printoptions(precision = 2)
np.set_printoptions(linewidth = 240)
#np.set_printoptions(suppress = True)

def jac_old(f, x):
	""" Computes the finite-difference approximation of the jacobian A=f'(x) of f at x"""
	fx = f(x)
	M = fx.size
	N = x.size
	A = sc.zeros((M,N))
	eps =  sqrt(finfo(double).eps)
	for i in range(0,N):
		x2 = copy(x)
		if x2[i] != 0:
			x2[i] = x2[i]*(1 + eps)
		else:
			x2[i] = eps
		fx2 = f(x2)
		A[:,i] = (fx2-fx)/eps
	return A

def jac(f, x, eps = sqrt(finfo(double).eps)):
	""" Computes the finite-difference approximation of the jacobian A=f'(x) of f at x"""
	fx = f(x)
	M = fx.size
	N = x.size
	A = sc.zeros((M,N))
	#eps =  sqrt(finfo(double).eps)
	for i in range(0,N):
		xp = copy(x)
		xm = copy(x)
		if x[i] != 0:
			xp[i] = xp[i]*(1 + eps)
			xm[i] = xm[i]*(1 - eps)
		else:
			xp[i] = eps
			xm[i] = -eps
		fxp = f(xp)
		fxm = f(xm)
		A[:,i] = (fxp-fxm)/(xp[i]-xm[i])
	return A

def f(u):
	return oderhs(u,0)[ind.xdot():]

def Sjoe(t,y):
	return oderhs(y,t)
Sjoe.name = 'Sjoe'

## Initial position/motion of the rolling elements
n = 3
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)

for i in xrange(0,n):
	u0[ind.r(i)] 		= (1+0*1e-8)*(b0+c0)/2.	# r
	u0[ind.theta(i)] 	= -pi/2.+2.*pi*i/n	# theta
	u0[ind.rdot(i)] 	= 0             	# rdot
	u0[ind.thetadot(i)] = b0/(b0+c0)*w0 	# thetadot
	u0[ind.phidot(i)] 	= b0/2./a0*w0    	# phidot - note that phi is not a state variable since it does not matter here


num = 16
J = zeros((num,11,19))
for i in range(0,num):
	J[i] = jac(f,u0,10**(-i))
for i in range(0,num): print i, max(max(J[i],0),0)
for k in range(0,num): print k, rank(J[k])

