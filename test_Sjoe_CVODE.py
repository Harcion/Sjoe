from Sjoe_m_func2 import *
from Sjoe_index import Index
from numpy import *
import IVP
import lsodar

def Sjoe(t,y):
	return oderhs(y,t)
Sjoe.name = 'Sjoe'

## Initial position/motion of the rolling elements
n = 8
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)

for i in xrange(0,n):
	u0[ind.r(i)] 		= (b0+c0)/2       # r
	u0[ind.theta(i)] 	= -pi/2 + 2*pi*i/n# theta
	#u0[ind.theta(i)] 	= -pi/2
	u0[ind.rdot(i)] 	= 0               # rdot
	u0[ind.thetadot(i)] = b0/(b0+c0)*w0   # thetadot
	u0[ind.phidot(i)] 	= b0/2/a0*w0      # phidot - note that phi is not a state variable since it does not matter here


u0 = array([ -1.76242809e-06,   2.76352696e-05,   2.49863406e-02,
		3.60926024e+01,   2.49894620e-02,   3.68814017e+01,
		2.49990368e-02,   3.76634995e+01,   2.50188206e-02,
		3.84516617e+01,   2.50287876e-02,   3.92372570e+01,
		2.50226201e-02,   4.00211008e+01,   2.50039592e-02,
		4.08050833e+01,   2.49913725e-02,   4.15901103e+01,
		-4.14597895e-03,  -1.02965855e-04,  -3.90691649e-04,
		3.77207067e+02,   1.88758319e+03,   1.95269428e-03,
		3.77162759e+02,   1.88700742e+03,   3.60118751e-03,
		3.77028471e+02,   1.88525136e+03,   4.92153159e-03,
		3.76198409e+02,   1.88837187e+03,   7.63066963e-04,
		3.75738352e+02,   1.89020148e+03,  -3.82863656e-03,
		3.75922411e+02,   1.88953549e+03,  -6.20071534e-03,
		3.76598617e+02,   1.88698681e+03,  -2.50644247e-03,
		3.77128621e+02,   1.88662007e+03])

cv=IVP.CVode(Sjoe,u0)
cv.set_method('BDF','Newton')
atol = array(zeros(N))
for i in range(0,N):
	atol[i] = 1e-8
cv.set_tolerance(atol,rtol=1.e-8)

numsteps = 1000
tend = 0.1
T = linspace(0,tend,numsteps)

cv(tend,numsteps)

stats=cv.stats(pr=1)
print stats

endstate = cv.aus[-1]

print u0
print endstate
print ""

for i in range(0,n):
	print endstate[ind.r(i)]
