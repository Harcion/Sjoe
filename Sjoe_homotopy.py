from math import *
from scipy import *
from pylab import *
from Sjoe_m_func import *

## Initial position/motion of the rolling elements
n = 8
N = 5*n+4
z1 = array(zeros(N))

for i in xrange(0,n):
	z1[5*i]   = (b0+c0)/2       # r
	z1[5*i+1] = 2*pi*i/n        # theta
	z1[5*i+2] = 0               # rdot
	z1[5*i+3] = b0/(b0+c0)*w0   # thetadot
	z1[5*i+4] = b0/2/a0*w0      # phidot - note that phi is not a state variable since it does not matter here


def H(u):
	z = array(zeros(N))
	
	for i in xrange(0,n):
		z[5*i]   = u[3*i]	# r
		z[5*i+3] = u[3*i+1]	# thetadot
		z[5*i+4] = u[3*i+2]	# phidot

	z[-2] = u[-2] # y
	
	l = u[-1] # Do something with this parameter
	
	t = 0 # Time?
	
	dz = oderhs(z, t)
	Hu = array(zeros(2+2*n))
	
	for i in xrange(0,n):
		Hu[2*i] = dz[5*i+3]		# thetadbldot
		Hu[2*i+1] = dz[5*i+4]	# phidbldot
	
	Hu[-2] = dz[-4]
	Hu[-1] = dx[-2]
		
u1 = array(zeros(3*n+2))