# %%
import matplotlib.pyplot as plt 
import numpy as np 

x = np.linspace(-10,10,100)
y = np.sin(x)
z = np.cos(x)
plt.plot(x,y)
plt.figure()
plt.plot(x,z)
plt.figure()
plt.plot(np.arange(100))
plt.show()




# %%
