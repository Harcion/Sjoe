from math import *
import numpy as np
import scipy as sc
from scipy import *
from pylab import *
from Sjoe_m_func import *

np.set_printoptions(precision = 2)
np.set_printoptions(linewidth = 180)
np.set_printoptions(suppress = True)

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
	z = u[:9]
	dz = oderhs(z,0)
	fu = array(zeros(5))
	fu[0:3] = dz[2:5]-u[5:8]
	fu[3:5] = dz[-2:]-dz[-4:-2]
	return fu
	


## Initial position/motion of the rolling elements
n = 1
N = 5*n+4
u0 = array(zeros(N))

for i in xrange(0,n):
	u0[5*i]   = (b0+c0)/2       # r
	u0[5*i+1] = 2*pi*i/n        # theta
	u0[5*i+2] = 0               # rdot
	u0[5*i+3] = b0/(b0+c0)*w0   # thetadot
	u0[5*i+4] = b0/2/a0*w0      # phidot - note that phi is not a state variable since it does not matter here

#print oderhs(u0,0)

#J = jac(lambda x: oderhs(x,0) , u0)
#
#print det(J)
#
#(q,r) = qr(J)
#(U, S, V) = svd(J)

#eps = 1e-8
#i = 1
#u1 = copy(u0)
#u1[i] = u1[i] + eps
#J1 = jac(lambda x: oderhs(x,0) , u1)
#print norm(J1-J)/eps


u0_ = array(zeros(14))
u0_[0:5] = u0[0:5]
u0_[5:8] = array(zeros(3))
u0_[8:12] = u0[-4:]
u0_[12:] = array(zeros(2))

print u0_

print f(u0_)

J = jac(f, u0_)
print J

(M,N) = shape(J)

(q,r) = sc.linalg.qr(J)
(U, s, Vh) = sc.linalg.svd(J)
S = sc.linalg.diagsvd(s,M,N)


print "U:", U

print "S:", S
print "V:", Vh

print norm(J - dot(U, dot(S, Vh)))
