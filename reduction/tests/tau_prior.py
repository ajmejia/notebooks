import scipy.stats as st
from scipy.optimize import root
from scipy.integrate import trapz as int_
import matplotlib.pyplot as plt
import numpy as np

x  = np.linspace(0, 2, 1000)
xx = np.linspace(0, 6, 1000)

pdf = lambda x, mu, sigma : st.norm.pdf(x, mu, sigma)

mu = 1.2

f = lambda sigma : int_(pdf(x, mu, sigma), x) - 0.68

result = root(f, 0.8)
sigma  = result.x[0]

print sigma
print result.success
print result.message
print int_(pdf(x, mu, sigma), x)
print int_(pdf(xx, mu, sigma), xx)

plt.plot(xx, pdf(xx, mu, sigma), "-", lw = 2, color = "#1A1A1A")
#l1, l2, = plt.plot(xx, pdf(xx, mu, sigma), "-", xx, pdf(xx, mu, 1.2724), "-", lw = 2)
#plt.legend([l1, l2], ["1st option", "2nd option"], frameon = False)
#plt.xlabel("tau")
#plt.ylabel("pdf")
#plt.tight_layout()
plt.show()
