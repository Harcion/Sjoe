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

np.set_printoptions(precision = 2)
np.set_printoptions(linewidth = 240)
#np.set_printoptions(suppress = True)

#def jac(f, x):
#	""" Computes the finite-difference approximation of the jacobian A=f'(x) of f at x"""
#	fx = f(x)
#	M = fx.size
#	N = x.size
#	A = sc.zeros((M,N))
#	eps =  sqrt(finfo(double).eps)
#	for i in range(0,N):
#		x2 = copy(x)
#		if x2[i] != 0:
#			x2[i] = x2[i]*(1 + eps)
#		else:
#			x2[i] = eps
#		fx2 = f(x2)
#		A[:,i] = (fx2-fx)/eps
#	return A

def jac(f, x):
	""" Computes the finite-difference approximation of the jacobian A=f'(x) of f at x"""
	fx = f(x)
	M = fx.size
	N = x.size
	A = sc.zeros((M,N))
	eps =  sqrt(finfo(double).eps)
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
		A[:,i] = (fxp-fxm)/(2.*eps)
	return A

#def f(u):
#	z = u[:ind.xdbldot()]
#	dz = oderhs(z,0)
#	fu = array(zeros(2+3*n))
#	fu = dz[ind.xdot():]-u[ind.xdbldot():]
#	return fu
	
def f(u):
	return oderhs(u,0)
	
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

#u0 = u0 + 1e-5*rand(N)

# Integrate first...
#
#cv=IVP.CVode(Sjoe,u0)
#cv.set_method('BDF','Newton')
#atol = array(zeros(N))
#for i in range(0,N):
#	atol[i] = 1e-8
#cv.set_tolerance(atol,rtol=1.e-8)
#
#numsteps = 1000
#tend = 0.03
#T = linspace(0,tend,numsteps)
#
#cv(tend,numsteps)
#
#stats=cv.stats(pr=1)
#print stats
#
u0_old = copy(u0)
#u0 = cv.aus[-1]  # Set new start state
u0 = array([    0.  ,     0.  ,     0.02,     9.72,     0.02,    11.82,     0.03,    13.92,    -0.02,     0.02,    -0.  ,   376.78,  1885.39,    -0.01,   377.34,  1888.48,     0.02,   375.34,  1894.02])
print u0

z0 = array(zeros(N-n-2))

z0[0] = 0 #x
z0[1] = 0 #y
k = 2
for i in xrange(0,n):
	z0[k] = u0[ind.r(i)]		; k+=1
	z0[k] = u0[ind.theta(i)]	; k+=1
	z0[k] = u0[ind.thetadot(i)]	; k+=1
	z0[k] = u0[ind.phidot(i)]	; k+=1
	
	
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
#
##print G(z0)

J = jac(G,z0)
(M,N) = J.shape

(Q,r,p) = qrp(J)
#print r
#print
# Set up the permutation matrix
P = zeros((N,N))
P[p,arange(N)] = 1

k = 2+3*n
R = r[:,:k]
S = r[:,k:]
V11 = -solve(R,S)

V = dot(P,r_[V11, eye(n)])  # Base for null space of J
print V
#print scipy.linalg.norm(dot(J,V))


#	
#z = copy(z0)
#print G(z)
#y = zeros(n)
#for i in range(0,n):
#	y[i] = u0[ind.phidot(i)]
#
#for i in range(0,10):
#	J = jac(G,z)
#	(M,N) = J.shape
#
#	(Q,r,p) = qrp(J)
#	#print p
#	#print
#	# Set up the permutation matrix
#	P = zeros((N,N))
#	P[p,arange(N)] = 1
#
#	k = 2+3*n
#	R = r[:,:k]
#	S = r[:,k:]
#	V11 = -solve(R,S)
#
#	V = dot(P,r_[V11, eye(n)])  # Base for null space of J
#	#print V
#	#print scipy.linalg.norm(dot(J,V))
#
#	zh = dot(V,y)
#	zp=- dot(pinv(J), G(z))
#
#	dz = zh+zp
#	#print scipy.linalg.norm(dot(J, dz) + G(z))
#	print 'J z_p',scipy.linalg.norm(dot(J, zp) + G(z),2)
#	print 'J pinv(J)',scipy.linalg.norm(dot(J, pinv(J)),2)
#	z = z + dz
#	#print G(z)
