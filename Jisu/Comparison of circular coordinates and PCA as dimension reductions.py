#!/usr/bin/env python
# coding: utf-8

# # Circular coordinates as a dimension reduction
# 
# The low-dimensional coordinate mappings by the circular
# coordinate is not only used as a visual representation, but also can
# be used as a data representation with reduced dimension so that it
# can be further feeded to data analysis frameworks. We first describe
# the pipeline for the circular coordinate as a dimension reduction,
# and compare it with PCA, in particular, how the topological structure
# of the data is preserved under these two dimension reductions. We
# will see with an example that linear dimension reduction methods such
# as the PCA can break down the topological structure of the data while
# the circular coordinate preserves the topological structure.
# 
# First, we describe how to use the circular coordinates
# as a dimension reduction. In other words, from the data $X$, we use
# the circular coordinates to create a data representation $X^{cc}$
# with reduced dimension so that $X^{cc}$ can be further feeded to
# data analysis frameworks. Given the circular coordinates $\theta_{1},\ldots,\theta_{k}:X=\{X_{1},\ldots,X_{n}\}\to S^{1}$,
# we understand $S^{1}\subset\mathbb{R}^{2}$ and create a map $\Theta:X\to\mathbb{R}^{2k}$
# as for each $X_{i}\in X$, 
# $$
# \Theta(X_{i})=(\theta_{1}(X_{i}),\ldots,\theta_{k}(X_{i}))\in\mathbb{R}^{2k}.
# $$
# Then the data representation $X^{cc}$ with reduced dimension is defined
# as $X^{cc}=\{\Theta(X_{1}),\ldots,\Theta(X_{n})\}$, i.e. $X_{i}^{cc}=\Theta(X_{i})$.
# 
# Using the circular coordinates as a dimension reduction
# has several benefits. First, since the reduced representation of the
# data is from the circular coordinates, which is again from the persistent
# cohomology, this representation focuses more on the topological structure
# of the data and less on how it is geometrically embedded. In particular,
# the circular coordinates representation does not depend on rotations
# or translations. Second, if we choose the circular coordinates $\theta_{1},\ldots,\theta_{k}$
# corresponding to significant topological features, the resulting circular
# coordinates representation focuses on these significant topological
# features and denoises topological noise in the data structure.

# In[1]:


get_ipython().magic(u'pylab inline')


# In[2]:


from Python_code import examples as eg
import numpy as np
from numpy import *
import dionysus


# # Dataset and its persistent cohomology
# Now, we compare how the topological structure of
# the data is preserved under the circular coordinates and the PCA.
# To see this, we use the dataset $X\subset\mathbb{R}^{3}$ in the plot,
# which is a uniform $60$ samples on a loop that wraps around the
# cube $[-1,1]^{3}$.

# In[3]:


from mpl_toolkits.mplot3d import Axes3D
m = 10
n = 12 * m

random.seed(0)

X1 = np.zeros((m, 3))
X1[:,0] = 1
X1[:,2] = np.array([uniform(0,1) for x in range(m)])

X2 = np.zeros((m, 3))
X2[:,2] = 1
X2[:,0] = np.array([uniform(0,1) for x in range(m)])

X3 = np.zeros((m, 3))
X3[:,2] = 1
X3[:,1] = np.array([uniform(0,1) for x in range(m)])

X4 = np.zeros((m, 3))
X4[:,1] = 1
X4[:,2] = np.array([uniform(0,1) for x in range(m)])

X5 = np.zeros((m, 3))
X5[:,1] = 1
X5[:,0] = - np.array([uniform(0,1) for x in range(m)])

X6 = np.zeros((m, 3))
X6[:,0] = -1
X6[:,1] = np.array([uniform(0,1) for x in range(m)])

X0 = np.vstack((X1, X2, X3, X4, X5, X6))
X = np.vstack((X0, -X0))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs = X[:,0], ys = X[:,1], zs = X[:,2])


# The persistent cohomology of $X$ is in the plot. As expected from the dataset, the $1$-dimensional persistent cohomology
# of $X$ has one prominent topological feature.

# In[4]:


prime = 23 #choose the prime base for the coefficient field that we use to construct the persistence cohomology.
vr = dionysus.fill_rips(X, 2, 4.) #Vietoris-Rips complex
cp = dionysus.cohomology_persistence(vr, prime, True) #Create the persistent cohomology based on the chosen parameters.
dgms = dionysus.init_diagrams(cp, vr) #Calculate the persistent diagram using␣ the designated coefficient field and complex.
dionysus.plot.plot_bars(dgms[1], show=True)
dionysus.plot.plot_diagram(dgms[1], show=True)
#dionysus.plot.plot_diagram(dgms[0], show=True)
#Plot the barcode and diagrams using matplotlib incarnation within Dionysus2. This mechanism is different in Dionysus.


# # Computing circular coordinates
# Then, we compute the circular coordinates from the longest cohomological feature in the persistent
# cohomology.

# In[5]:


threshold = 0.5
bars = [bar for bar in dgms[1] if bar.death-bar.birth > threshold] #choosing cocycle that persist at least threshold=0.5
cocycles = [cp.cocycle(bar.data) for bar in bars]
chosen_cocycle= cocycles[0]
chosen_bar= bars[0]


# In[6]:


vr_8 = dionysus.Filtration([s for s in vr if s.data <= max([bar.birth for bar in bars])])
coords = dionysus.smooth(vr_8, chosen_cocycle, prime)


# In[7]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs = X[:,0], ys = X[:,1], zs = X[:,2], c = coords)


# # Circular coordinates representation and its persistent cohomology
# For the circular coordinates representation, we let
# $\theta:X\to S^{1}\subset\mathbb{R}^{2}$ to be the circular coordinate
# from the prominent cobundary of the persistent cohomology in the previous plot,
# and create a map $\Theta:X\to\mathbb{R}^{2}$ by this single circular
# coordinate, i.e. $\Theta(X_{i})=\theta(X_{i})$. The resulting circular
# coordinates representation $X^{cc}=\{\Theta(X_{1}),\ldots,\Theta(X_{n})\}$
# is in the plot.

# In[8]:


X_cc = numpy.zeros((size(coords), 2))
X_cc[:, 0] = cos(2 * pi * numpy.matrix(coords))
X_cc[:, 1] = sin(2 * pi * numpy.matrix(coords))

scatter(X_cc[:, 0], X_cc[:, 1])
plt.axis('equal')
plt.title('Circular Coordinates')
plt.show()


# And the persistent cohomology of $X^{cc}$ is in the plot. As you can see, the 1-dimensional persistent cohomology of $X^{cc}$
# in the plot contains one topological feature, which is the same with
# the 1-dimensional persistent cohomology of the original data $X$.
# Hence the topological structure of the original data is
# preserved.

# In[9]:


prime = 23 #choose the prime base for the coefficient field that we use to construct the persistence cohomology.
vr_cc = dionysus.fill_rips(X_cc, 2, 4.) #Vietoris-Rips complex
cp_cc = dionysus.cohomology_persistence(vr_cc, prime, True) #Create the persistent cohomology based on the chosen parameters.
dgms_cc = dionysus.init_diagrams(cp_cc, vr_cc) #Calculate the persistent diagram using␣ the designated coefficient field and complex.
dionysus.plot.plot_bars(dgms_cc[1], show=True)
dionysus.plot.plot_diagram(dgms_cc[1], show=True)
#dionysus.plot.plot_diagram(dgms[0], show=True)
#Plot the barcode and diagrams using matplotlib incarnation within Dionysus2. This mechanism is different in Dionysus.


# # PCA representation and its persistent cohomology
# For the PCA representation, to be comparable with
# the circular coordinates, we choose the projection dimension to be
# $2$. The resulting PCA representation $X^{pca}$ is in the plot. From
# the plot, we can see that the $1$-dimensional cohomological structure
# of the original data $X$ is distorted in $X^{pca}$.

# In[10]:


from Python_code import PCAtool as PCAtool
X_pca = PCAtool.pca(X, K=2)

scatter(X_pca[:, 0], X_pca[:, 1])
plt.axis('equal')
plt.title('PCA')
plt.show()


# For the PCA representation, to be comparable with
# the circular coordinates, we choose the projection dimension to be
# $2$. The resulting PCA representation $X^{pca}$ is in the plot. From
# the plot, we can see that the $1$-dimensional cohomological structure
# of the original data $X$ is collapsed in $X^{pca}$. And the persistent
# cohomology of $X^{pca}$ is in the plot. As you can see, the 1-dimensional
# persistent cohomology of $X^{pca}$ in the plot contains $3$ topological features,
# which is different from the 1-dimensional persistent cohomology of
# the original data $X$. Hence the topological structure
# of the original data is not preserved.

# In[11]:


prime = 23 #choose the prime base for the coefficient field that we use to construct the persistence cohomology.
vr_pca = dionysus.fill_rips(X_pca, 2, 4.) #Vietoris-Rips complex
cp_pca = dionysus.cohomology_persistence(vr_pca, prime, True) #Create the persistent cohomology based on the chosen parameters.
dgms_pca = dionysus.init_diagrams(cp_pca, vr_pca) #Calculate the persistent diagram using␣ the designated coefficient field and complex.
dionysus.plot.plot_bars(dgms_pca[1], show=True)
dionysus.plot.plot_diagram(dgms_pca[1], show=True)
#dionysus.plot.plot_diagram(dgms[0], show=True)
#Plot the barcode and diagrams using matplotlib incarnation within Dionysus2. This mechanism is different in Dionysus.

