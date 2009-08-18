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

np.set_printoptions(precision = 2)
np.set_printoptions(linewidth = 240)
#np.set_printoptions(suppress = True)

def jac(f, x):
	""" Computes the finite-difference approximation of the jacobian A=f'(x) of f at x"""
	fx = f(x)
	M = fx.size
	N = x.size
	A = sc.zeros((M,N))
	eps = 1e-8
	for i in range(0,N):
		e = sc.zeros(N)
		e[i] = eps
		x2 = x + e
		fx2 = f(x2)
		A[:,i] = (fx2-fx)/eps
	return A

def f(u):
	z = u[:ind.xdbldot()]
	dz = oderhs(z,0)
	fu = array(zeros(2+3*n))
	fu = dz[ind.xdot():]-u[ind.xdbldot():]
	return fu
	


## Initial position/motion of the rolling elements
n = 3
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)

for i in xrange(0,n):
	u0[ind.r(i)] 		= (1+1e-8)*(b0+c0)/2.	# r
	u0[ind.theta(i)] 	= -pi/2.+2.*pi*i/n	# theta
	u0[ind.rdot(i)] 	= 0             	# rdot
	u0[ind.thetadot(i)] = b0/(b0+c0)*w0 	# thetadot
	u0[ind.phidot(i)] 	= b0/2./a0*w0    	# phidot - note that phi is not a state variable since it does not matter here

z0 = array(zeros(N-n-2))

z0[0] = 0 #x
z0[1] = 0 #y
k = 2
for i in xrange(0,n):
	z0[k] = u0[ind.r(i)] ;	k+=1
	z0[k] = u0[ind.theta(i)] ;	k+=1
	z0[k] = u0[ind.thetadot(i)] ;	k+=1
	z0[k] = u0[ind.phidot(i)] ;	k+=1
	
	
def G(z):
	n = (z.size-2)/4
	u = array(zeros(5*n+4))
	u[ind.x()] = z[0]
	u[ind.y()] = z[1]
	for i in xrange(0,n):
		u[ind.r(i)] 		= z[2+4*i]	# r
		u[ind.theta(i)] 	= z[3+4*i]	# theta
		u[ind.rdot(i)] 	    = 0         # rdot
		u[ind.thetadot(i)]  = z[4+4*i] 	# thetadot
		u[ind.phidot(i)] 	= z[5+4*i]  # phidot - note that phi is not a state variable since it does not matter here
	
	du = oderhs(u,0)
	return du[ind.xdot():]

#print G(z0)

J = jac(G,z0)
(M,N) = J.shape

k = 2+3*n

(Q,r,p) = qrp(J)

# Set up the permutation matrix
P = zeros((N,N))
P[p,arange(N)] = 1

R = r[:,:k]
S = r[:,k:]
V11 = -solve(R,S)

V = dot(P,r_[V11, eye(n)])

print V
#print scipy.linalg.norm(dot(J,V))
