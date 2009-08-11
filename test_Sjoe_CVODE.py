from Sjoe_m_func2 import *
from Sjoe_index import Index
from numpy import *
import IVP
import lsodar

def Sjoe(t,y):
	return oderhs(y,t)
Sjoe.name = 'Sjoe'
## Initial position/motion of the rolling elements
n = 1
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)

for i in xrange(0,n):
	u0[ind.r(i)] 		= (b0+c0)/2       # r
	u0[ind.theta(i)] 	= 2*pi*i/n        # theta
	u0[ind.rdot(i)] 	= 0               # rdot
	u0[ind.thetadot(i)] = b0/(b0+c0)*w0   # thetadot
	u0[ind.phidot(i)] 	= b0/2/a0*w0      # phidot - note that phi is not a state variable since it does not matter here


cv=IVP.CVode(Sjoe,u0)
cv.set_method('BDF','Newton')
atol = array(zeros(N))
for i in range(0,N):
	atol[i] = 1e-8
cv.set_tolerance(atol,rtol=1.e-8)

numsteps = 1000
tend = 0.03
T = linspace(0,tend,numsteps)

cv(tend,numsteps)

stats=cv.stats(pr=1)
print stats

print u0
print cv.aus[-1]

(Y, info) = lsodar.odeintr(oderhs, copy(u0), T, atol = 1e-8, rtol = 1e-8, full_output=1)

print Y[-1,:]

print cv.aus[-1] - Y[-1,:]
#	del cv
#	assert stats['Steps']==577
