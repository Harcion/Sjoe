from math import *
import numpy as np
import scipy as sc
import pylab as pylab
from scipy import *
from pylab import *
from Sjoe_m_func2 import *
from Sjoe_index import *

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
	z = u[:ind.xdbldot()]
	dz = oderhs(z,0)
	fu = array(zeros(2+3*n))
	fu = dz[ind.xdot():]-u[ind.xdbldot():]
	return fu
	


## Initial position/motion of the rolling elements
n = 4
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)

for i in xrange(0,n):
	u0[ind.r(i)] 		= (b0+c0)/2			# r
	u0[ind.theta(i)] 	= -pi/2+2*pi*i/n	# theta
	u0[ind.rdot(i)] 	= 0             	# rdot
	u0[ind.thetadot(i)] = b0/(b0+c0)*w0 	# thetadot
	u0[ind.phidot(i)] 	= b0/2/a0*w0    	# phidot - note that phi is not a state variable since it does not matter here

u0_ = r_[u0, array(zeros(2+3*n))]
print u0_
print f(u0_)
print ""

J = jac(f, u0_)
#print J

(M,N) = shape(J)

Js = array(zeros(J.shape))
#exc = array([ind.xdot(), ind.ydot(), ind.xdbldot(), ind.ydbldot()])
#exc = array([ind.x(), ind.y(), ind.xdot(), ind.ydot()])
exc = array([ind.xdot(), ind.ydot(), ind.xdbldot(), ind.ydbldot()])
q = array([ind.x(), ind.y()])

for i in range(0,n):
	#exc = r_[exc, array([ind.rdot(i), ind.thetadot(i), ind.rdbldot(i), ind.thetadbldot(i), ind.phidbldot(i)])]
	#exc = r_[exc, array([ind.r(i), ind.theta(i), ind.rdot(i), ind.thetadot(i), ind.phidot(i)])]
	exc = r_[exc, array([ind.rdot(i), ind.thetadbldot(i), ind.phidbldot(i)])]
	q = r_[q, array([ind.r(i), ind.theta(i), ind.thetadot(i), ind.phidot(i), ind.rdbldot(i)])]
	
k = 0
#for i in set(range(0,N)) - set(exc):
#	Js[:,k] = copy(J[:,i])
#	k += 1
#
#for i in exc:
#	Js[:,k] = copy(J[:,i])
#	k += 1

for i in q:
	Js[:,k] = copy(J[:,i])
	k += 1

for i in set(range(0,N)) - set(q):
	Js[:,k] = copy(J[:,i])
	k += 1



#P  = Js[:,:2+3*n]
P  = Js[:,:q.size]  #(N-exc.size)
print P.shape
#print P
prank = pylab.rank(P)
print "Rank:", prank
max_rank = 2+3*n
print "Max rank:", max_rank
print "Js rank:", pylab.rank(Js)

#(MP, NP) = P.shape
#(U, s, Vh) = sc.linalg.svd(P)
#S = sc.linalg.diagsvd(s, MP, NP)
#
#print U
#print ""
#print S
#print ""
#print Vh
#
#print pylab.rank(S)

#
#for i in range(0,N):
#	P2 = c_[P , Js[:,i]]
#	if pylab.rank(P2) > prank:
#		print i
#		for k in range(0,N):
#			P3 = c_[P2 , Js[:,k]]
#			if pylab.rank(P3) == max_rank:
#				print "   ", k


#for i in range(0,max_rank):
#	P3 = copy(P)
#	P3 = delete(P3, i, axis=1) # delete ith column
#	if pylab.rank(P3) == prank:
#		print i
		
#Js[:, ind.phidot(0)],  Js[:, ind.rdbldot(3)], Js[:, ind.rdbldot(1)],
#for i in range(0,N):
#	P3 = c_[P2 , Js[:,i]]
#	if pylab.rank(P3) == max_rank:
#		print i
