import scipy.stats as st
from scipy.optimize import root
from scipy.integrate import trapz as int_
import matplotlib.pyplot as plt
import numpy as np

x  = np.linspace(0.1, 1.0, 1000)
xx = np.linspace(0.0, 1.0, 1000)

pdf = lambda x, mu, sigma : st.norm.pdf(x, mu, sigma)

mu = 0.3

f = lambda sigma : int_(pdf(x, mu, sigma), x) - 0.68

result = root(f, 0.31)
sigma  = result.x[0]

print sigma
print result.success
print result.message
print int_(pdf(x, mu, sigma), x)
print int_(pdf(xx, mu, sigma), xx)

plt.plot(xx, pdf(xx, mu, sigma), "-", lw = 2, color = "#4D4D4D")
plt.show()

# Ivan's realization of mu pdf ---------------------------------------------------------------------
#import random as ran
#
#def draw_any(dom,dist,comp):
	#if len(dom) == len(dist) and len(dom) == len(comp):
		#not_success = True
		#while not_success:
			#i = ran.choice(xrange(len(dist)))
			#u = np.random.uniform(0,1,1)
			#if u < dist[i]/comp[i]:
				#not_success = False
				#return dom[i]
	#else:
		#print 'dom, dist and comp must have the same length'
		#sys.exit(1)
#
#mu_grid = 1000
#dist_mean = 0.3
#dist_sd = 1.272367 #<-- ivan
#dist_sd = 0.36566883 #<-- alfredo
#
#domain = np.linspace(-10,10,mu_grid)
#
#dist_mu = (1/(dist_sd*np.sqrt(2*np.pi))) * np.exp((-(domain-dist_mean)**2)/(2*dist_sd**2)) #Prior PDF shape
#A_bad = np.sum(dist_mu)
#dist_mu = dist_mu/A_bad #Prior PDF normalization
#
#comp_dist = np.ones(mu_grid)*dist_mu[np.argmax(dist_mu)]*1.1 #Comparison distribution --Caution with factor--

