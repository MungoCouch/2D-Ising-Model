import numpy as np 
from matplotlib import pyplot as plt#
from itertools import product

	
class Ising:
	def __init__(self,s,T):

		self.size = (s,s) # tuple (height,width) of lattice.
		self.Lattice = np.random.choice([-1,1],self.size)
		self.temperature = T
		self.sites = list(product(range(s),range(s)))
	#def magnetization():
	def RelativeProbability(self,DelE): 
		
		return np.exp(-DelE/self.temperature)
	#def TotalEnergy(): return None
	def mag(self): return sum(self.Lattice[s[0]][s[1]]for s in self.sites)/(self.size[0]**2)
	def LocalEnergy(self,i,j):
		return -self.Lattice[i][j]*(self.Lattice[(i+1)%self.size[0]][j] + self.Lattice[i][(j+1)%self.size[0]] + self.Lattice[(i-1)%self.size[0]][j] + self.Lattice[i][(j-1)%self.size[0]])
	def totalEnergy(self):
		return sum([self.LocalEnergy(i,j) for (i,j) in self.sites])
	def NewState(self): 
		([i],[j]) = np.random.randint(self.size[0],size = (2,1))
		DelE = -self.LocalEnergy(i,j)
		Pr = self.RelativeProbability(DelE)
		return (i,j),Pr,DelE
	# Given M_N samples
	# N num of samples
	# k number of samples to test correlation.
	def autoCor(self,M,N,k):
		""" Given M_N samples
	 	N num of samples
	 	k number of samples to test correlation.
	 	we compute the corresponding auto-correlation.
		"""
		E_M = (1/N)*sum(M)	
		Var_M = (1/N)*sum([(m_i - E_M)**2 for m_i in M])
		Num = (1/(N-k))*sum([(M[i] - E_M)*(M[k+i]) for i in range(N-k)])
		return Num/Var_M
	def Tau(self,M,N):

		"""autocorrelation time"""
		return sum([(1-(k/N))*self.autoCor(M,N,k) for k in range(N-1)])
	def Met_1(self):

		(i,j),Pr,DelE = self.NewState() # generate a new possible state
		if DelE < 0: 
			self.Lattice[i,j] *= -1
			return None
		#generate a uniform random number to compare against acceptance ratio
		#If the new state is favourable again the random number change state
		#otherwise discard new state
		u = np.random.uniform(0,1) 
		#print(Pr,u)
		if min([1,Pr]) >= u:
			
			self.Lattice[i,j] *= -1
		return None

	def Met(self, EachSweep = None):
		"""For the purposes of analysis we may want 
		to take data points more regularly than once 
		per sweeep, but the default is 1 per sweep of the whole lattice"""
		if EachSweep:
			EachSweep = int(EachSweep) 
			for i in range(EachSweep): 
				self.Met_1()
			return None 
		#If only one spin site is updated per MC then one cannot hope for uncorrelated samples
		for i in range(self.size[0]**2): 
			self.Met_1()
		return None 
