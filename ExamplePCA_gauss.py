from Python_code import examples as eg
import numpy as np
from numpy import *
import dionysus
from Python_code import PCAtool as PCAtool
def denoise_dgm(dgms,threshold=0.1):
    threshold = .1
    bars_0 = [np.array(bar).tolist() for bar in dgms[0] if bar.death-bar.birth > threshold] 		#choosing cocycle that persist at least threshold=1.
    bars_1 = [np.array(bar).tolist() for bar in dgms[1] if bar.death-bar.birth > threshold] 		#choosing cocycle that persist at least threshold=1.
    #cocycles = [cp.cocycle(bar.data) for bar in bars]
    #plt is the matplotlib incarnation.
    #print( type(bars_1) )
    #bars_1=np.array(bars_1)
    #print(type(bars_1[1]))
    #print(bars_1[1])
    #print(type(bars_1))

    dgm_denoise_0=bars_0
    dgm_denoise_1=bars_1
    #print(type(dgm_denoise_0))
    #dionysus.plot.plot_diagram(dgm_denoise_0, show=False)
    #dionysus.plot.plot_diagram(dgm_denoise_1, show=False)
    return([dgm_denoise_0,dgm_denoise_1])

#Let us compute the bottle-neck distance between full dataset and k-PCA dataset. k denotes the number of principal components.
prime=23
D=20
#Here is the data insertion part
mean=np.repeat(0,D)
cov=np.identity(D)
dat1=np.random.multivariate_normal(mean, cov, 100)
vr = dionysus.fill_rips(dat1, 2, 10) #Vietoris-Rips complex
cp = dionysus.cohomology_persistence(vr, prime, True) #Create the persistent cohomology
dgms = dionysus.init_diagrams(cp, vr) #Calculate the persistent diagram using the designated coefficient field.
#record for plottings..
A=np.array([0,0,1,1])

for thres in [0,0.025,0.05,0.1,0.2]:
    #######Thresholding the features on the persistent diagrams describing persistent cohomology of VR complex(CC)
    dgms_tmp=denoise_dgm(dgms,thres)#thresholding CC diagrams
    dgms_base_0=dionysus._dionysus.Diagram()#take persistent points of dim 0
    dgms_base_1=dionysus._dionysus.Diagram()#take persistent points of dim 1
    #pybind issue, need to be correct instead of succint.
    for itr in range(len(dgms_tmp[0])):
        dgms_base_0.append(dgms_tmp[0][itr])
    for itr in range(len(dgms_tmp[1])):
        dgms_base_1.append(dgms_tmp[1][itr])
    #print(type(dgms_base[0]))
    #dionysus.bottleneck_distance(dgms_base_0, dgms_base_1)
    #######
    for k in range(D):
        dat_k = PCAtool.pca(dat1,K=k+1)
        vr_k = dionysus.fill_rips(dat_k, 2, 10) #Vietoris-Rips complex
        cp_k = dionysus.cohomology_persistence(vr_k, prime, True) #Create the persistent cohomology
        dgms_k = dionysus.init_diagrams(cp_k, vr_k) #Calculate the persistent diagram using the designated coefficient field.
        #######Thresholding the features on the persistent diagrams describing persistent cohomology of kPCA-VR complex(kPCA)
        #############Here we use the same thresholding for fairness. But this configuration is easy to change.
        PCA_thres=thres
        #############
        dgms_k_tmp=denoise_dgm(dgms_k,PCA_thres)#thresholding CC diagrams
        dgms_k_0=dionysus._dionysus.Diagram()#take persistent points of dim 0
        dgms_k_1=dionysus._dionysus.Diagram()#take persistent points of dim 1
        #pybind issue, need to be correct instead of succint.
        for itr in range(len(dgms_k_tmp[0])):
            dgms_k_0.append(dgms_k_tmp[0][itr])
        for itr in range(len(dgms_k_tmp[1])):
            dgms_k_1.append(dgms_k_tmp[1][itr])
        #print(type(dgms_base[0]))
        #dionysus.bottleneck_distance(dgms_base_0, dgms_base_1)
        #######
        bdist_0 = dionysus.bottleneck_distance(dgms_base_0, dgms_k_0)
        bdist_1 = dionysus.bottleneck_distance(dgms_base_1, dgms_k_1)
        #dionysus.plot.plot_diagram(dgms_k[1], show=True) 
        newrow=[k+1,thres, bdist_0,bdist_1]
        print(newrow)
        A = np.vstack([A, newrow])
print(A)
#Column names of A
#1 number of principal components k; 
#2 threshold used for filtering out topological features by persistence;
#3 Bottleneck distance between dimension 0 diagrams of CC and PCA
#3 Bottleneck distance between dimension 1 diagrams of CC and PCA
np.savetxt('ExamplePCA_gauss_100.txt',A)
