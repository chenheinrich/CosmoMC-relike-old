import numpy as np
import GetTau
from matplotlib import pyplot as plt

import getdist.plots as gplot
g = gplot.getSinglePlotter(chain_dir=r'/project/kicp/chenhe/Reionization/cosmomc/chains/test6_chains/')
  
samples = g.sampleAnalyser.samplesForRoot('test6')

p = samples.getParams()

# process parameters
z1 = 0.0
z2 = 30.0
nz = 5

reion_nbasis = 5
reion_mj = np.empty(reion_nbasis)

#npt = p.H0.size
npt = 3
cumtau_all = np.empty((npt, nz))
tau_30 = np.empty(npt)

for ip in np.arange(0,npt): 
	print ' '
	print 'ip = ', ip, ', npt = ', npt
	H0 = p.H0[ip]
	h = H0/100.0
	h2 = h**2
	omegab = p.omegabh2[ip]/h2
	omegac = p.omegach2[ip]/h2
	omegam = p.omegam[ip]
	omegal = p.omegal[ip]
	#omegak = p.omegak[ip]
	omegak = 0.0
	#yheused = p.yheused[ip]
	yheused = 0.245341
	reion_mj[0] = p.reion_m1[ip]
	reion_mj[1] = p.reion_m2[ip]
	reion_mj[2] = p.reion_m3[ip]
	reion_mj[3] = p.reion_m4[ip]
	reion_mj[4] = p.reion_m5[ip]
	#reion_m1 = -0.113946226566
	#reion_m2 = -0.137704508707
	#reion_m3 = 0.261631502869
	#reion_m4 = -0.239478806731
	#reion_m5 = 0.069993393982
	#reion_mj[0] = reion_m1
	#reion_mj[1] = reion_m2
	#reion_mj[2] = reion_m3
	#reion_mj[3] = reion_m4
	#reion_mj[4] = reion_m5
	
	# make sure that Basis.zmax_basis == 30.0
	tau_30[ip] = GetTau.GetTauMax(H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj)
	print 'tau_30[ip] = ', tau_30[ip]
	
	[z2_array, cumtau] = GetTau.GetTauSeq(z1, z2, nz, H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj)
	cumtau_all[ip,:] = cumtau

print 'tau_30 = ', tau_30
print 'cumtau_all = ', cumtau_all

#for iz in np.arange(0, nz):
#	derived = cumtau_all[:,[iz]]
#	mean[iz] = samples.mean(derived)
#	err[iz] = samples.std(derived)
#	limit95[iz] = samples.twoTailLimits(derived, 0.95)
#	samples.addDerived(derived, name='tau%f'%nz[iz], label=r'\tau_{%f}'%nz[iz])
	
#samples.updateBaseStatistics()

#print mean
#print err
#print limit95

#fig, ax = plt.subplots(1,figsize=(6,4.5))
#plt.errorbar(z2_array, mean, yerr=err, fmt='.')

#legloc1 = 'upper right'
#title1 = r'$\tau\mathrm{\ vs\ }z_{max}$'
#xlabel1 = r'$z_{max}$'
#ylabel1 = r'$\tau(0, z_{max})$'

#plt.legend(loc=legloc1)	
#plt.title(title1)
#ax = plt.gca()
#ax.set_xlabel(xlabel1)	
#ax.set_ylabel(ylabel1)	
#ax.set_xlim(np.array([z1, z2]))


	