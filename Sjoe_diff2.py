from math import *
import numpy as np
import scipy as sc
from scipy import *
from pylab import *
from Sjoe_m_func2 import *

np.set_printoptions(precision = 2)
np.set_printoptions(linewidth = 180)
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
	z = u[:-5*n]
	dz = oderhs(z,0)
	fu = array(zeros(5*n))
	fu = dz[-5*n:]-u[-5*n:]
	return fu
	


## Initial position/motion of the rolling elements
n = 1
N = 5*n+4
u0 = array(zeros(N))

for i in xrange(0,n):
	u0[2+2*i] 		= (b0+c0)/2       # r
	u0[3+2*i] 		= 2*pi*i/n        # theta
	u0[4+2*n + 3*i] = 0               # rdot
	u0[5+2*n + 3*i] = b0/(b0+c0)*w0   # thetadot
	u0[6+2*n + 3*i] = b0/2/a0*w0      # phidot - note that phi is not a state variable since it does not matter here

u0_ = r_[u0, array(zeros(5*n))]

#print u0_
#print f(u0_)

J = jac(f, u0_)
print J

(M,N) = shape(J)

#(q,r) = sc.linalg.qr(J)
#(U, s, Vh) = sc.linalg.svd(J)
#S = sc.linalg.diagsvd(s,M,N)
#
#print "U:", U
#print "S:", S
#print "V:", Vh
#print sc.linalg.norm(J - dot(U, dot(S, Vh)))

Js = array(zeros(J.shape))
exc = array([2+2*n, 3+2*n, 4+5*n, 5+5*n])
# xdot, ydot, xdoubledot, ydoubledot
#exc = array([0,1, 2+2*n, 3+2*n])
# x, y, xdot, ydot
for i in range(0,n):
	#exc = r_[exc, array([4+2*n+3*i, 5+2*n+3*i, 6+2*n+3*i, 4+5*n + 3 + 5*i, 4+5*n + 4 + 5*i])]#, 2+2*i, 3+2*i])]
	# rdot, thetadot, phidot, thetadoubledot, phidoubledot  # r, theta
	
	#exc = r_[exc, array([2+2*i, 3+2*i, 4+2*n+3*i, 5+2*n+3*i, 6+2*n+3*i])]#, ])]
	# r, theta, rdot, thetadot, phidot
	
	#exc = r_[exc, array([4+2*n+3*i, 4+5*n + 2 + 5*i, 4+5*n + 3 + 5*i, 4+5*n + 4 + 5*i])]#, 2+2*i, 3+2*i])]
	# rdot, thetadot, phidot, rdoubledot, thetadoubledot, phidoubledot  # r, theta

	exc = r_[exc, array([4+2*n+3*i, 5+2*n+3*i, 4+5*n + 2 + 5*i, 4+5*n + 3 + 5*i, 4+5*n + 4 + 5*i])]#, 2+2*i, 3+2*i])]
	# rdot, thetadot, phidot, rdoubledot, thetadoubledot, phidoubledot  # r, theta

k = 0
for i in set(range(0,N)) - set(exc):
	Js[:,k] = copy(J[:,i])
	k += 1

for i in exc:
	Js[:,k] = copy(J[:,i])
	k += 1

print ""
#print Js

print Js[:5*n,:5*n]
print rank(Js[:5*n,:5*n])
#
#def permute(P, frm, to):
#	(M,N) = P.shape
#	if (frm == to) or (frm == N + to) or (to == frm + N):  # Don't need to permute anything
#		print frm
#		print to
#		return P
#	P[frm, to]  = 1
#	P[to, frm]  = 1
#	P[frm,frm] = 0
#	P[to,to] = 0
#	return P
#
#P = diag(ones(N)) # Permutation matrix
#P = permute(P, 2+2*n, -1) # xdot
#P = permute(P, 3+2*n, -2) # ydot
#P = permute(P, 4+5*n, -3) # xdoubledot
#P = permute(P, 5+5*n, -4) # ydoubledot
#P = permute(P, 4+2*n, -5) # rdot
#P = permute(P, 4+5*n + 2, -6) # thetadoubledot
#P = permute(P, 4+5*n + 3, -7) # phidoubledot