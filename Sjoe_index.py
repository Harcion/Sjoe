class Index:
	class BodyOutOfBounds(Exception):
		def __init__(self, value = 'Body does not exist'):
			self.value = value
		def __str__(self):
			return repr(self.value)
	
	def __init__(self, num_rollers):
		self.n = num_rollers

	def x(self):
		return 0
	
	def y(self):
		return 1
	
	def r(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 2 + 2*i
	
	def theta(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 3 + 2*i
	
	def xdot(self):
		return 2 + 2*self.n
	
	def ydot(self):
		return 3 + 2*self.n
	
	def rdot(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 4 + 2*self.n + 3*i
	
	def thetadot(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 5 + 2*self.n + 3*i
	
	def phidot(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()		
		return 6 + 2*self.n + 3*i
	
	def xdbldot(self):
		return 4 + 5*self.n
	
	def ydbldot(self):
		return 5 + 5*self.n
	
	def rdbldot(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 6 + 5*self.n + 3*i
	
	
	def thetadbldot(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 7 + 5*self.n + 3*i
	
	def phidbldot(self, i):
		if i < 0 or i >= self.n:
			raise self.BodyOutOfBounds()
		return 8 + 5*self.n + 3*i	