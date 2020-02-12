import scipy as sp
import numpy as np

def pca(X,K=1,centering=1):
  '''
  Transform a data matrix of size N by M into its first K principal coordinates.

  Parameters
  ----------
  X: float() - N by M matrix, each row is a datapoint lying in R^M
  K: int() - number of points
  centering: boolean() - whether impose an auto-centering procedure before PCA.

  Output
  ------
  np.array of dimensions (N,M)
  '''
  N, M = X.shape
  if centering==1:
      X_bar=np.mean(X, axis = 0)
      X=X-X_bar
  if K>M:
      print('PCA Error: Number of principal components cannot exceed the dimension of data.')
  Sigma = np.dot(X.T, X) / (N-1)
  eigen_vals, eigen_vecs = np.linalg.eig(Sigma)
  #Retain the eigen_vecs corresponding to the K largest eigen_vals.
  eigen_ord = np.argsort(eigen_vals)
  #Extract the indices ordered according to descending magnitude.
  eigen_ord = eigen_ord[0:K]
  #Retain only largest K.
  # Project X onto PC space
  X_pca = np.dot(X, eigen_vecs[:,eigen_ord])
  return X_pca#eigen_ord, eigen_vecs, eigen_vals
