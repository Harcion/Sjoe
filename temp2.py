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
from scipy.linalg import norm

np.set_printoptions(precision = 2)
np.set_printoptions(linewidth = 240)
#np.set_printoptions(suppress = True)

def jac_old(f, x, eps =  sqrt(finfo(double).eps)):
	""" Computes the finite-difference approximation of the jacobian A=f'(x) of f at x"""
	fx = f(x)
	M = fx.size
	N = x.size
	A = sc.zeros((M,N))
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
	for i in range(0,N):
		xp = copy(x)
		xm = copy(x)
		
		if x[i] != 0:
			h = abs(x[i]*eps)
		else:
			h = eps
			
		xp[i] += h
		xm[i] -= h
		
		fxp = f(xp)
		fxm = f(xm)

		A[:,i] = (fxp-fxm)/(2*h)#/sqrt(h)
		
	return A

def f(u):
	return oderhs(u,0)#[ind.xdot():]
	
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

eps = 1e-3
max_it = 1000
max_half = 10
h = 1

acc0 = f(u0)
acc0_ = acc0[ind.xdot():]
u = copy(u0)
u += acc0/norm(acc0)

acc1 = f(u)
acc1_ = acc1[ind.xdot():]
r = zeros(max_it)
err = zeros(max_it)

i = 0
while i < max_it and norm(acc1_) > eps:
	d0 = acc0_[0:2]/norm(acc0_)
	
	same_direction = dot(d0,acc1_[0:2])
	#print "%e" % u[ind.r(0)]
	if  same_direction < 0:
		h = h/2.
		u = copy(u0)
		u += h*acc0/norm(acc0)
		acc1 = f(u)
		acc1_ = acc1[ind.xdot():]
		print "Halving,", h, i
		print acc1
	else:
		u0 = u
		u = copy(u) + h*acc1/norm(acc1)
		acc0 = acc1
		acc0_ = acc0[ind.xdot():]
		acc1 = f(u)
		acc1_ = acc1[ind.xdot():]
		#print acc1
		
#		if h < 1:
#			h *= 1.1
	
		i += 1
	#print f(u)[ind.ydot()]
