import numpy as np
from matplotlib import pyplot as plt

#import getdist.plots as gplot
# change this to just reading in the chains myself, not relying on GetDist.
#g = gplot.getSinglePlotter(chain_dir=r'/project/kicp/chenhe/Reionization/cosmomc/chains/')
#samples = g.sampleAnalyser.samplesForRoot('test6')
#p = samples.getParams()
#m1 = p.reion_m1
#m2 = p.reion_m2
#m3 = p.reion_m3
#m4 = p.reion_m4
#m5 = p.reion_m5
#liketot = (p.chi2_lowTEB + p.chi2_plik + p.chi2_prior)/2.0

chain_dir = '/project/kicp/chenhe/Reionization/cosmomc/chains/test6_chains/'
root = 'test6'
i = 1
sample_tmp = np.genfromtxt(chain_dir + root + '_%s'%i + '.txt')
for i in np.arange(2,5):
	sample_tmp = np.genfromtxt(chain_dir + root + '_%s'%i + '.txt')
	sample = np.vstack((sample, sample_tmp))
	
# mjs: cols (5 to 9) + 2, index (4 to 8) + 2
sample_like = sample[:,1] 
indp_m1 = 4
reion_nbasis = 5
indc_m1 = indp_m1+2
m1 = sample[:,indc_m1 ]  # to delete (for convenience right now
m2 = sample[:,indc_m1 +1] 	
#print liketot[0:5]
print sample_like[0:5]

# import covariance matrix
in_dir = 'results/test6/'
covmat = np.genfromtxt(in_dir + root + '.covmat')
cov = covmat[indp_m1:indp_m1+reion_nbasis,indp_m1:indp_m1+reion_nbasis]	
print cov

#ngrid_1d = 100
#grid_m1 = np.linspace(np.min(m1), np.max(m1), ngrid_1d)
#grid_m2 = np.linspace(np.min(m1), np.max(m1), ngrid_1d)
#m1_sparse, m2_sparse = np.meshgrid(grid_m1, grid_m2)
#z = np.sin(m1_sparse**2 + m2_sparse**2)
#h = plt.contourf(grid_m1, grid_m2, z)
frac_sigma = 1.0
sample_weight = sample[:,0] 
ms_sample = sample[:,indc_m1:indc_m1+reion_nbasis]

i = 0
ngrid = np.empty(reion_nbasis)
for ind in np.arange(indc_m1, indc_m1+reion_nbasis):
	#width = frac_sigma * np.std(sample[:,ind])
	width = frac_sigma * cov[i,i]
	diff = np.max(sample[:,ind]) - np.min(sample[:,ind])
	ngrid[i] = diff/width
	i = i + 1
	
ngrid = np.floor(ngrid)
print ngrid	

H, edges = np.histogramdd(ms_sample, bins = (ngrid), weights = sample_weight, normed = True)
grid_m1 = edges[0]
grid_m2 = edges[1]



#g = gplot.getSinglePlotter(chain_dir=chain_dir)
#def my_loglike_fun(m1, m2):
#	return np.sin(m1**2 + m2**2)
#loglike = my_loglike_func(m1_sparse, m2_sparse)
#density = Density2D(grid_m1,grid_m2, np.exp(-loglike / 2))
#density.contours = np.exp(-np.array([1.509, 2.4477]) ** 2 / 2)
#g.add_2d_contours(root, 'reion_m1', 'reion_m2', filled=True, density=density)
	
