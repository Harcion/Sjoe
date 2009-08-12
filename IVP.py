# -*- coding: UTF-8 -*-

from __future__ import division
#import pylab, numpy, scipy, scipy.linalg
import numpy as np
import scipy as sp
import pylab as PL
from numpy import array, dot, linspace, zeros, ones, shape, arange, log10, exp, sqrt, pi, hstack, column_stack, vstack, floor, vectorize, ndarray
from numpy.linalg import solve, norm, eig
from pylab import plot, xlabel, ylabel, grid, axis, legend, xlabel, ylabel, title, subplot, imshow, colorbar, text
from scipy import meshgrid
from scipy.linalg.matfuncs import expm
from scipy.interpolate import interp1d


import sundials_cv as cv

import nose.tools as nt
import numpy.testing as npt



class IVPSolver(object):
	"""
	Generic class for initial value problem solvers
	Olivier Verdier, Claus Führer
	Lund University
	"""


	def __init__(self, f, u0, t0=0):
		self.f = f   # rhs function
		self.u0= u0  # initial value
		self.t0= t0  # initial time
		self.ts = [t0]
		self.us = [u0]


	def adjust_stepsize(self):
		pass
	
	# default values 
	h = .01  # initial step size (not used in CVODE)
	tf = 1.      # final time
	max_steps = 10000 # max no of steps 
	rtol=1.e-4   # relative tolerance
	atol=1.e-4   # absolute tolerance
	
	def generator(self, t, u, tf,nt):
		"""
		Generates the (t,u) values until t > tf
		"""
		for i in xrange(self.max_steps):
			if t >= tf:
				break
			t, u = self.step(t, u)
			yield t,u
			self.adjust_stepsize()
			self.h=min(self.h,abs(tf-t))
		else:
			raise Exception("Final time not reached within max_steps steps")
		
	
	def __call__(self, tfinal,nt=None):
		# start from the last time we stopped
		t = t0 = self.ts[-1]
		u = self.us[-1]
		tus = list(self.generator(t, u, tfinal,nt))
		
		self.ts.extend(q[0] for q in tus)
		self.us.extend(q[1] for q in tus)
		
		self.ats = array(self.ts)
		self.aus = array(self.us)

	def plot(self):
		"""
		Plot the computed solution.
		"""
		if not hasattr(self,'ats'):
			raise Exception, 'No data to plot.'
		plot(self.ats, self.aus, '.-')
		xlabel('time')
		ylabel('state')
		try:
			title(self.f.name)
		except:
			pass

	def plot2D(self):
		"""
		Plot ux vs uy
		"""
		plot(self.aus[:,0],self.aus[:,1], '.-')
		xlabel('ux')
		ylabel('uy')
	
	quiver_res = 20
	def quiver(self):
		mins = self.aus.min(axis=0)
		maxs = self.aus.max(axis=0)
		X,Y = np.meshgrid(linspace(mins[0], maxs[0], self.quiver_res), 
								linspace(mins[1], maxs[1], self.quiver_res))
		Z = np.dstack([X,Y])
		vals = self.f(0,Z.transpose(2,0,1))
		PL.quiver(X,Y,vals[0], vals[1])
	
	def compute_error(self):
		"""Compare u with the exact solution at current time"""
		t = self.ts[-1]
		u = self.us[-1]
		return norm(u - self.f.exact(self.us[0], self.ts[0], t)) 

	def error_for_h(self, h):
		"""
		Compute the global error for a given stepsize h.
		"""
		self.h = h
		self.run()
		error = self.compute_error()
		return error
	
	error_resolution = 10
	err_kmin=2;err_kmax=3
	def plot_error(self, do_plot=True):
		"""
		Plot the total error with respect to the number of steps.
		
		Parameters
		----------
		do_plot : boolean(True)
			do plot the error points (set to False mostly for testing purposes)
		
		Class Parameters
		----------------
		err_kmin, err_kmax : number
			bounds of the exponents of N to be used.
		error_resolution : number
			number of points to use.
		
		Returns
		-------
		slope : number
			the slope of the linear regression
		self.reg : tuple
			the result of the regression is stored for possible later use
		"""
		k_logs = linspace(self.err_kmin, self.err_kmax, self.error_resolution)
		hs = self.time/(10**k_logs)
		error_logs = log10(array([self.error_for_h(h) for h in hs]))
		self.reg = linear_regression(k_logs, error_logs, do_plot)
		slope = self.reg[0]
		if do_plot:
			axis('equal')
			plot(k_logs, error_logs, 'o', label='computed error')
			title('log error vs. log N (%s)' % type(self).__name__)
			xlabel("log10(N)")
			ylabel("log10(error)")
			legend()
			grid(True)
		return slope
	
	
	def plot_steps(self):
		ats = self.ats
		hs = np.diff(ats)
		plot((ats[1:] + ats[:-1])/2, hs)

class ExplicitEuler (IVPSolver):
	def step(self, t, u):
		return t + self.h, u + self.h*self.f(t, u)


class RungeKutta4 (IVPSolver):
	"""
	Runge-Kutta of order 4.
	"""
	def step(self, t, u):
		f = self.f
		h = self.h
		Y1 = f(t, u)
		Y2 = f(t + h/2., u + h*Y1/2.)
		Y3 = f(t + h/2., u + h*Y2/2.)
		Y4 = f(t + h, u + h*Y3)
		return t+h, u + h/6.*(Y1 + 2.*Y2 + 2.*Y3 + Y4)

class RungeKutta34 (IVPSolver):
	"""
	Adaptive Runge-Kutta of order four.
	Obs. Step rejection not implemented.
	"""
	_order = 4.
	# default tolerance

	def adjust_stepsize(self):
		self.h *= (self.atol/self.error)**(1/self._order)

	def step(self, t, u):
		f = self.f
		h = self.h
		Y1 = f(t, u)
		Y2 = f(t + h/2., u + h*Y1/2.)
		Y3 = f(t + h/2, u + h*Y2/2)
		Z3 = f(t + h, u - h*Y1 + 2*h*Y2)
		Y4 = f(t + h, u + h*Y3)
		self.error = norm(h/6*(2*Y2 + Z3 - 2*Y3 - Y4))
		return t+h, u + h/6*(Y1 + 2*Y2 + 2*Y3 + Y4)

class CVode(IVPSolver):
	"""
	Class CVODE   Initial value problem solver based on BDF or Adams
	multistep method
	initializes and calls the serial, dense version of CVODE from the 
	SUNDIALS (2.4.0) package.
	
	(read the SUNDIALS/CVODE documentation for method details
	and in particular for default values.)
	
	Claus Führer, Numerical Analysis, Lund University
	"""
	
	def __init__(self, f, u0, t0=0):
		super(CVode, self).__init__(f,u0,t0)
		self.CVdata=cv.Cvode_wrap(len(u0))
		# defaults
		self.CVdata.abstol=self.atol
		self.CVdata.reltol=self.rtol
		self.CVdata.abstol_ar=np.zeros(self.CVdata.dim)
		self.CVdata.max_steps=self.max_steps
		self.maxord=None
	def set_method(self,discr='ADAMS',iter='Fixed Point'):
		if discr=='BDF':
			self.CVdata.discr=2
		if iter=='Newton':
			self.CVdata.iter=2
	def set_tolerance(self,atol,rtol):
		if not isinstance(atol,float):
			self.CVdata.abstol_ar=array(atol)
		else:
			self.CVdata.abstol_ar=np.zeros(self.CVdata.dim)
			self.abstol=atol
		if not isinstance(rtol,float):
			raise Exception,"CVODE relative tolerance must be a (scalar) float."
		self.CVdata.reltol=rtol
	def set_max_ord(self,maxord):
		"""
		Sets the maximal order of the method:
		defaults = maximal values:
		Adams:  maxord=12
		BDF  :  maxord= 5
		
		An input value greater than the default will result in the default value.
		"""
		self.maxord=maxord
	def eval_rhs(self,t,u):
		return cv.eval_rhs(t,u,self.f)
	def generator(self,t,u,tfinal,nt=100):
		self.CVdata.cvinit(t,self.f,u,self.maxord)
		return self.CVdata.run(t,tfinal,nt)
	def stats(self,pr=None):
		stats=self.CVdata.stats
		if pr==1:
			try:
				problem=self.f.name
			except NameError:
				problem='------'
			print 'Integration statistics for Problem %s \n\n' % problem
			print stats
		return stats
	

	




