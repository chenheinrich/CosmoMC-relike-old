import numpy as np
import ReionTau

def GetTauMax(H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj):
	
	Tau = ReionTau.SetupTau(H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj)
	tau = Tau(0.0, 30.0) # need to change this if Basis.zmax_basis changes

	return tau
	
def GetTau(z1, z2, H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj):
	
	Tau = ReionTau.SetupTau(H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj)
	tau = Tau(z1, z2)

	return tau	

def GetTauSeq(z1, z2, nz, H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj):
	
	nz = np.floor(nz) # round down nz
	if ((z1 <= z2) & (nz >= 1)):
		dz = (z2 - z1) / nz
		z_array = np.arange(z1, z2+dz, dz)
		z1_array = z_array[:-1]
		z2_array = z_array[1:]
			
		Tau = ReionTau.SetupTau(H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj)
		tau = np.empty(z1_array.size)
		i = 0
		for z1, z2 in zip(z1_array, z2_array):	
			tau[i] = Tau(z1, z2)
			i = i + 1
		cumtau = np.cumsum(tau)
	else:
		print 'Error from GetTauSeq: z1 > z2 or/and nz < 1'		
	return z2_array, cumtau	