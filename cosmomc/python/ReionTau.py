import ReionXe
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate

# ------- Setup once for all points on the chain ------- 

# physical constants (mass_ratio_He_H also set in ReionXe, make sure consistent)
# cosmology settings
# files: neutrino evolution
# keys: USE_XEZ_CAMB (also needed in ReionXe), tau_zmin_basis = 0.03875651

# Physical constants (output from camb in SI unit, cross checked with http://pdg.lbl.gov/1998/consrpp.pdf)
sigma_thomson   =   0.665246160000E-28
kappa     		=   0.167730888412E-08 #8piG
c         		=   0.299792458000E+09
sigma_boltz   	=   0.567040000000E-07
m_H       		=   0.167357500000E-26
Mpc		  		=   0.308567800000E+23
mass_ratio_He_H =   0.397150000000E+01

# Physical parameters 
tcmb	 		=   0.272550000000E+01
nu_massless_degeneracy = 3.046 * 2.0/3.0
# 1 massive neutrino by default

# File names 
fn_nu_evol = 'python/tanh_base2p6_grhonu.dat' # temporary, might change

# omeganuh2 evolution in z output from camb: z = 3 to 40
nu_evol = np.genfromtxt(fn_nu_evol) # 2 cols: z vs grhonu_mass(z) (used in this code)
print 'nu_evol has z range: %f to %f'%(nu_evol[0,0], nu_evol[-1,0])		
#omeganu_evol = nu_evol[:,1]*(1.0+nu_evol[:,0])**4/3.0/(H0**2.0)
#omeganu_evol_interpolate = interpolate.interp1d(nu_evol[:,0], omega_nu_evol, kind='cubic')
grhonu_mass_interpolate = interpolate.interp1d(nu_evol[:,0], nu_evol[:,1], kind='cubic')
grhonu_mass_eval = np.vectorize(grhonu_mass_interpolate)

# integration related
zmin_basis = ReionXe.Basis.zmin_basis
zmax_basis = ReionXe.Basis.zmax_basis
tau_zmin_basis = 0.03875651 # tau(0 to 6.0) for pc	
USE_XEZ_CAMB = False


# ------- Setup once for all points on the chain ------- 

def GetGrho(H0, omegab, omegac, omegam, omegal, omegak): # # grho = 8piG * rho/c^2
	grhom = 3.0*(H0/c*1000.0)**2.0 #3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
	grhog = kappa/c**2.0*4.0*sigma_boltz/c**3.0*tcmb**4.0*Mpc**2.0 #8*pi*G/c^2*4*sigma_B/c^3 T^4
	# grhog=1.4952d-13*tcmb**4
	grhor = 7.0/8.0*(4.0/11.0)**(4.0/3.0)*grhog #7/8*(4/11)^(4/3)*grhog (per neutrino species)
	# grhor=3.3957d-14*tcmb**4
	grhornomass = grhor*nu_massless_degeneracy   
	grhoc = grhom*omegac
	grhob = grhom*omegab
	grhol = grhom*omegal
	grhok = grhom*omegak
	return grhom, grhog, grhob, grhoc, grhornomass, grhol, grhok 	
				
class SetupTau():
	
	def __init__(self, H0, omegab, omegac, omegam, omegal, omegak, yheused, reion_nbasis, reion_mj):
		self.omegab = omegab
		self.yheused = yheused
		[self.grhom, self.grhog, self.grhob, self.grhoc, self.grhornomass, self.grhol, self.grhok] = GetGrho(H0, omegab, omegac, omegam, omegal, omegak)
		self.yheused = yheused
		self.reion_nbasis = reion_nbasis
		self.reion_mj = reion_mj
		self.limit_nstep = 1000
		
	def GetDtda(self, z):
		a = 1.0/(1.0 + z)
		a2 = a**2
		#  8*pi*G*rho*a**4.
		grhoa2 = self.grhok*a2 + (self.grhoc + self.grhob) *a + self.grhog + self.grhornomass + self.grhol*a2**2
		grhoa2=grhoa2 + grhonu_mass_eval(z)
		dtda = np.sqrt(3.0/grhoa2)
		return dtda
		
	def GetAkthom(self):
		Nnow = self.omegab*(1.0-self.yheused)*self.grhom*c**2.0/kappa/m_H/Mpc**2.0
		akthom = sigma_thomson * Nnow * Mpc
		return akthom
		
	def Integrand(self, z):	
		if USE_XEZ_CAMB:
			integ = -ReionXe.xe_camb_eval(z) * self.GetDtda(z)
			#integ = -xe_camb_eval(z) * dtauda_eval(z)
		else:
			integ = ReionXe.xe_pc_eval(z, self.reion_nbasis, self.reion_mj) * self.GetDtda(z) # pcs
			#integ = xe_tanh_eval(z) * dtda(z) # tanh
		return integ
	
	def __call__(self, z1, z2):
		#print 'zmin_basis = ', zmin_basis, ', zmax_basis = ', zmax_basis
		print ''
		warn_tol = 1e-5
		if (z2 <= zmax_basis):
			if (z1 == 0.0): # expect this case mostly
				print 'Integrating tau from: zmin_basis = ', zmin_basis, ' to z2 = ', z2
				int = integrate.quad(self.Integrand, zmin_basis, z2, limit=self.limit_nstep)
				akthom = self.GetAkthom()
				tau = akthom * int[0] + tau_zmin_basis
				if (int[0] != 0):
					if (int[1]/int[0] > warn_tol):
						print 'Warning: relative error of integral exceeds %e.'%warn_tol
				#print int
				print 'added to tau_zmin_basis = ', tau_zmin_basis
			elif (z1 <= zmax_basis): # just integrate, special case (maybe for testing)
				print 'Integrating tau from: z1 = ', z1, ' to z2 = ', z2
				int = integrate.quad(self.Integrand, z1, z2, limit=self.limit_nstep)
				if (int[0] != 0):
					if (int[1]/int[0] > warn_tol):
						print 'Warning: relative error of integral exceeds %e.'%warn_tol
				#print int
				akthom = self.GetAkthom()
				tau = akthom * int[0]
			else:
				print 'Error from SetupTau: integration lower limit z1 exceeds zmax_basis'
		else:
			print 'Error from SetupTau: integration upper limit z2 exceeds zmax_basis'			
		return tau
		
# dtauda from camb for testing: log in z from 4 to -13
#fn_dtauda = 'test8_5pc_dtauda.dat'
#dtauda = np.genfromtxt(fn_dtauda) # 2 cols: z vs xe(z) (used in this code)
#print 'dtauda has z range: %f to %f'%(dtauda[0,0],dtauda[-1,0])			
#dtauda_interpolate = interpolate.interp1d(dtauda[:,0], dtauda[:,1], kind='cubic')
#dtauda_eval = np.vectorize(dtauda_interpolate)
		
			