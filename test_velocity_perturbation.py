from Sjoe_m_func2 import *
from Sjoe_index import Index
from numpy import *
import IVP
import time

def Sjoe(t,y):
	return oderhs(y,t)
Sjoe.name = 'Sjoe'

## Initial position/motion of the rolling elements
n = 8
N = 5*n+4
u0 = array(zeros(N))
ind = Index(n)

#for i in xrange(0,n):
#	u0[ind.r(i)] 		= (b0+c0)/2       # r
#	u0[ind.theta(i)] 	= -pi/2 + 2*pi*i/n# theta
#	u0[ind.rdot(i)] 	= 0               # rdot
#	u0[ind.thetadot(i)] = b0/(b0+c0)*w0   # thetadot
#	u0[ind.phidot(i)] 	= b0/2/a0*w0      # phidot - note that phi is not a state variable since it does not matter here

#cv=IVP.CVode(Sjoe,u0)
#cv.set_method('BDF','Newton')
#atol = array(zeros(N))
#for i in range(0,N):
#	atol[i] = 1e-8
#cv.set_tolerance(atol,rtol=1.e-8)
#
#numsteps = 100
#t0 = 0.005
#
#ex_time0 = time.clock()
#cv(t0,numsteps)
#ex_time0 = time.clock() - ex_time0
#
#stats=cv.stats(pr=1)
#print stats
#
#Y = cv.aus
#t = cv.ts
#state0 = Y[-1]
#
#for i in range(0,n):
#	subplot(n+1,5,5*i+1) # r
#	plot(t, Y[:,ind.r(i)])
#
#	subplot(n+1,5,5*i+2) # th
#	plot(t, Y[:,ind.theta(i)])
#
#	subplot(n+1,5,5*i+3) # rdot
#	plot(t, Y[:,ind.rdot(i)])
#
#	subplot(n+1,5,5*i+4) # thetadot
#	plot(t, Y[:,ind.thetadot(i)])
#
#	subplot(n+1,5,5*i+5) # phidot
#	plot(t, Y[:,ind.phidot(i)])
#
#subplot(n+1,5,5*n+1)
#plot(t, Y[:,ind.x()])
#subplot(n+1,5,5*n+2)
#plot(t, Y[:,ind.xdot()])
#subplot(n+1,5,5*n+3)
#plot(t, Y[:,ind.y()])
#subplot(n+1,5,5*n+4)
#plot(t, Y[:,ind.ydot()])
#
#show()
#
#len0 = len(cv.ts)
#
#t1 = 2*t0
#
#ex_time1 = time.clock()
#cv(t1,numsteps)
#ex_time1 = time.clock() - ex_time1
#
#Y = cv.aus[len0:]
#t = cv.ts[len0:]
#state0 = Y[-1]
#
#figure()
#for i in range(0,n):
#	subplot(n+1,5,5*i+1) # r
#	plot(t, Y[:,ind.r(i)])
#
#	subplot(n+1,5,5*i+2) # th
#	plot(t, Y[:,ind.theta(i)])
#
#	subplot(n+1,5,5*i+3) # rdot
#	plot(t, Y[:,ind.rdot(i)])
#
#	subplot(n+1,5,5*i+4) # thetadot
#	plot(t, Y[:,ind.thetadot(i)])
#
#	subplot(n+1,5,5*i+5) # phidot
#	plot(t, Y[:,ind.phidot(i)])
#
#subplot(n+1,5,5*n+1)
#plot(t, Y[:,ind.x()])
#subplot(n+1,5,5*n+2)
#plot(t, Y[:,ind.xdot()])
#subplot(n+1,5,5*n+3)
#plot(t, Y[:,ind.y()])
#subplot(n+1,5,5*n+4)
#plot(t, Y[:,ind.ydot()])
#
#show()
#



#
## Endstate at t=0.005:
# array([ -3.12389648e-07,   2.82812046e-05,   2.50096813e-02,
#		3.14817787e-01,   2.50262817e-02,   1.09934532e+00,
#		2.50282524e-02,   1.88060072e+00,   2.50144220e-02,
#		2.66663069e+00,   2.49960013e-02,   3.45423783e+00,
#		2.49875915e-02,   4.24115807e+00,   2.49866347e-02,
#		5.02732275e+00,   2.49936367e-02,   5.81286402e+00,
#		-7.46800018e-03,  -7.59069491e-04,   2.77472966e-03,
#		3.76622261e+02,   1.88670729e+03,   8.29568845e-04,
#		3.76254001e+02,   1.88824350e+03,  -1.57856688e-03,
#		3.75725818e+02,   1.89094483e+03,  -3.08763092e-03,
#		3.76585387e+02,   1.88703174e+03,  -1.43424948e-03,
#		3.77056381e+02,   1.88573443e+03,  -4.33668989e-04,
#		3.77137066e+02,   1.88709100e+03,   8.07733445e-04,
#		3.77137365e+02,   1.88720102e+03,   1.57545715e-03,
#		3.77056246e+02,   1.88600149e+03])
#
## Endstate at t=0.01:
#array([ -7.71179565e-07,   2.77200897e-05,   2.50241704e-02,
#		2.19440825e+00,   2.50064179e-02,   2.98048553e+00,
#		2.49923984e-02,   3.76239529e+00,   2.49864966e-02,
#		4.55196442e+00,   2.49887012e-02,   5.34042404e+00,
#		2.49977697e-02,   6.12728920e+00,   2.50169088e-02,
#		6.91226752e+00,   2.50284747e-02,   7.69513877e+00,
#		3.00859059e-04,  -4.27336583e-04,  -6.39188912e-03,
#		3.75812567e+02,   1.88934890e+03,  -1.06498480e-02,
#		3.76828366e+02,   1.88567604e+03,  -4.37114190e-03,
#		3.77172556e+02,   1.88672359e+03,  -7.89437877e-04,
#		3.77293844e+02,   1.88799973e+03,   3.23063768e-03,
#		3.77240661e+02,   1.88748347e+03,   5.64121419e-03,
#		3.77047727e+02,   1.88549098e+03,   8.61246496e-03,
#		3.76267256e+02,   1.88751593e+03,   1.55778140e-03,
#		3.75710508e+02,   1.88946348e+03])
#
# Execution times:
# 0.005:  187.48000000000013
# 0.01 :   45.1400000000001



u0 = array([ -3.12389648e-07,   2.82812046e-05,   2.50096813e-02,
		3.14817787e-01,   2.50262817e-02,   1.09934532e+00,
		2.50282524e-02,   1.88060072e+00,   2.50144220e-02,
		2.66663069e+00,   2.49960013e-02,   3.45423783e+00,
		2.49875915e-02,   4.24115807e+00,   2.49866347e-02,
		5.02732275e+00,   2.49936367e-02,   5.81286402e+00,
		-7.46800018e-03,  -7.59069491e-04,   2.77472966e-03,
		3.76622261e+02,   1.88670729e+03,   8.29568845e-04,
		3.76254001e+02,   1.88824350e+03,  -1.57856688e-03,
		3.75725818e+02,   1.89094483e+03,  -3.08763092e-03,
		3.76585387e+02,   1.88703174e+03,  -1.43424948e-03,
		3.77056381e+02,   1.88573443e+03,  -4.33668989e-04,
		3.77137066e+02,   1.88709100e+03,   8.07733445e-04,
		3.77137365e+02,   1.88720102e+03,   1.57545715e-03,
		3.77056246e+02,   1.88600149e+03])



cv=IVP.CVode(Sjoe,u0)
cv.set_method('BDF','Newton')
atol = array(zeros(N))
for i in range(0,N):
	atol[i] = 1e-8
cv.set_tolerance(atol,rtol=1.e-8)

numsteps = 100

# Perturb the velocities:
for i in range(0,n):
	u0[ind.thetadot(i)] *= (1 + rand()/10.)


t1 = 0.005

ex_time1 = time.clock()
cv(t1,numsteps)
ex_time1 = time.clock() - ex_time1

Y = cv.aus
t = cv.ts
state0 = Y[-1]

figure()
for i in range(0,n):
	subplot(n+1,5,5*i+1) # r
	plot(t, Y[:,ind.r(i)])

	subplot(n+1,5,5*i+2) # th
	plot(t, Y[:,ind.theta(i)])

	subplot(n+1,5,5*i+3) # rdot
	plot(t, Y[:,ind.rdot(i)])

	subplot(n+1,5,5*i+4) # thetadot
	plot(t, Y[:,ind.thetadot(i)])

	subplot(n+1,5,5*i+5) # phidot
	plot(t, Y[:,ind.phidot(i)])

subplot(n+1,5,5*n+1)
plot(t, Y[:,ind.x()])
subplot(n+1,5,5*n+2)
plot(t, Y[:,ind.xdot()])
subplot(n+1,5,5*n+3)
plot(t, Y[:,ind.y()])
subplot(n+1,5,5*n+4)
plot(t, Y[:,ind.ydot()])

show()


##
# extime1 = 63.180000000000064
#           56.940000000000055
#			57.240000000000009
