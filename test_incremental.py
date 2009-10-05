from Sjoe_m_func2 import *
from Sjoe_index import Index
from numpy import *
import scipy
import numpy
import IVP
import time

def jac(f, x, eps = sqrt(finfo(double).eps)):
	""" Computes the finite-difference approximation of the Jacobian A=f'(x) of f at x"""
	fx = f(x)
	M = fx.size
	N = x.size
	A = numpy.zeros((M,N))
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

		A[:,i] = (fxp-fxm)/(2*h)

	return A

def Sjoe(args = None):
	if args != None:
		def rhs(t,y):
			return oderhs(t,y, args)
	else:
		def rhs(t,y):
			return oderhs(t,y)
	rhs.name = 'rhs'
	return rhs

## Initial position/motion of the rolling elements
n = 8
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)
# This is rand(8)*1e-6:
perturbation = array([  9.78702856e-07,   2.93684668e-07,   4.79000309e-07,
						3.21827339e-07,   8.01445882e-07,   7.26611584e-08,
						8.74822963e-08,   7.41679693e-07])

for i in xrange(0,n):
	u0[ind.r(i)]        = (b0+c0)/2 + perturbation[i] #rand()*1e-6       # r
	u0[ind.theta(i)]    = -pi/2 + 2*pi*i/n# theta
	u0[ind.rdot(i)]     = 0               # rdot
	u0[ind.thetadot(i)] = 0*b0/(b0+c0)#*w0   # thetadot
	u0[ind.phidot(i)]   = 0*b0/2/a0*628.32      # phidot - note that phi is not a state variable since it does not matter here

end_state = array([9.81, 628.32, 1000.])
h = 1. / numpy.max(end_state)

#rhs_ = Sjoe((0,0,0))
#rhs = Sjoe(10*h*end_state)
rhs_ = Sjoe((0.0981,0,0))
def f(u):
	return rhs_(0,u)
A0 = jac(f,u0)

uk = copy(u0)
#
#rhs_ = Sjoe((0.0981,0,0))
#def f(u):
#	return rhs_(0,u)
#
#f2 = f(uk)
#print f2-f1

#print A0

x_solve = numpy.linalg.solve(A0, f(u0))

x_lstsq = numpy.linalg.lstsq(A0, f(u0))[0]

res_solve = dot(A0, x_solve)-f(u0)

res_lstsq = dot(A0, x_lstsq)-f(u0)

print "norm(res_solve)", norm(res_solve)

print "norm(res_lstsq)", norm(res_lstsq)


print "scipy.linalg.svdvals(A0)", scipy.linalg.svdvals(A0)

print "norm(x_lstsq)", norm(x_lstsq)
print "norm(x_solve)", norm(x_solve)

print "numpy.max(x_lstsq)", numpy.max(x_lstsq)
print "numpy.max(x_solve)", numpy.max(x_solve)






#while(norm(f(uk)) > 5):
#	uk -= lstsq(jac(f,uk), f(uk))[0]
#
#print uk
#print f(uk)
#print norm(f(uk))
#A = jac(f,uk)
#
#rhs = Sjoe(h*end_state)
#print h*end_state
#def f(u):
#	return rhs(0,u)
#
#A02 = jac(f,u0)
#uk2 = copy(uk)
#
##for i in range(0,10):
##	uk2_norm = norm(f(uk2))
##	print uk2_norm
##	if uk2_norm < 10:
##		break
##	uk2 -= lstsq(jac(f,uk2), f(uk2))[0]
##
##print uk2
##print f(uk2)
