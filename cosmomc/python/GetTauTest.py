import scipy.integrate as integrate
import scipy.interpolate as interpolate
import numpy as np
from matplotlib import pyplot as plt 

# user selection 
zmin = 6.0
zmax = 30.0 # can set this based on basis

#zmax = 33.9102500352015 # 5pc
#zmin = 3.247064042977144E-006 # 5pc 

#zmax = 0.14810020247539E+02 # tanh
#zmin = 0.32470640429771E-05 # tanh

tau_zmax = -0.000142555165612 # tau(30 to 33.9..) for pc
tau_zmin = 0.0387633897876 # tau(3e-6 to 6.0) for pc
#tau_zmin = 0.03875651 # tau(0 to 6.0) for pc
z1 = zmin # user can select (put this in a function)
z2 = zmax # user can select
xe_fid = 0.15
print 'zmin, zmax are: ', zmin, zmax

# file names 
fn_xez_camb = 'test8_test_tau_integrator_xe/test8_5pc_xez_with_helium_corrected.dat' # cosmomc xe(z) used for testing
fn_xez_recomb = 'test8_test_tau_integrator_xe/test8_5pc_xez_recomb_long.dat' # cosmomc xe(z) recomb used for testing
fn_nu_evol = 'test8_test_tau_integrator_xe/test8_5pc_grhonu.dat'
fn_basis = 'xepcs.dat' # not original mortonson file
mortonson_format = True
INTERP_BASIS_ORDER = 1
USE_XEZ_CMB = True
USE_XEZ_RECOMB = True
LOAD_XEZ_RECOMB = True
INTERP_KIND_XEZ_CAMB = 'cubic'

# load from chain, different for each point
omegabh2 = 0.022338
omegach2 = 0.118685
omeganu2 = 0.000645
H0 = 0.677539062500E+02
omegam = 0.308494
omegal = 0.691506
yheused = 0.2453
omegak = 0.0
reion_nbasis = 5
reion_m1 = -0.113946226566
reion_m2 = -0.137704508707
reion_m3 = 0.261631502869
reion_m4 = -0.239478806731
reion_m5 = 0.069993393982

# process parameters
h = H0/100.0
h2 = h**2
omegab = omegabh2/h2
omegac = omegach2/h2
reion_mj = np.empty(reion_nbasis)
reion_mj[0] = reion_m1
reion_mj[1] = reion_m2
reion_mj[2] = reion_m3
reion_mj[3] = reion_m4
reion_mj[4] = reion_m5
#for j in np.arange(0, reion_nbasis):
#	reion_mj[j] = reion_m1

# constants from camb (SI unit; cross checked with http://pdg.lbl.gov/1998/consrpp.pdf)
sigma_thomson   =   0.665246160000E-28
kappa     		=   0.167730888412E-08 #8piG
c         		=   0.299792458000E+09
sigma_boltz   	=   0.567040000000E-07
m_H       		=   0.167357500000E-26
Mpc		  		=   0.308567800000E+23
mass_ratio_He_H =   0.397150000000E+01

# some cmb parameters that could be changed
tcmb	 		=   0.272550000000E+01
nu_massless_degeneracy = 3.046 * 2.0/3.0
is_cosmological_constant = True
num_nu_massive  = 1

#---------------- Read in basis & add points --------------

if mortonson_format:
	# assume original Mortonson file format
	# row = 1: z; row > 1: Bj(z)
	basis = np.genfromtxt(fn_basis) 
	basis = np.transpose(basis) 
	# col = 1: z; col > 1: Bj(z)
	# so jth basis is accessible with basis(:,j)
	zs = basis[:,0]
	nz = basis.shape[0] # number of rows
	nbasis = basis.shape[1] - 1 # number of columns - 1
	if zs[0] > zs[1]:
		print 'Error: need the z in %s to be monotonically increasing.'%fn_basis
else:
	# assume transposed Mortonson file format with one additional row
	# col = 1: z; col > 1: Bj(z)
	# so jth basis is accessible with basis(:,j)
	basis = np.genfromtxt(fn_basis, skip_header = 1) 
	zs = basis[0,:]
	nz = basis.shape[1] # number of cols
	nbasis = basis.shape[0] - 1 # number of rows - 1

# Extend grid by making two additional rows
zmin_basis = zs[0] - (zs[1] - zs[0]) # assuming increasing z order ! CAN CHECK
zmax_basis  = zs[-1] + (zs[-1] - zs[-2])
if (zs[1] - zs[0]) == (zs[-1] - zs[-2]):
	dz_basis = zs[1] - zs[0]
rowmin = np.hstack((zmin_basis , np.zeros(nbasis)))
rowmax = np.hstack((zmax_basis , np.zeros(nbasis)))	
# remake basis with added rows
basis = np.vstack((rowmin, basis, rowmax))
zs = basis[:,0]
nz = basis.shape[0]
nbasis = basis.shape[1] - 1 # number of columns - 1


print 'Done reading in and extending basis.'

#------------------ Setup basis -------------------------------  

basis_list = []
for j in np.arange(0, reion_nbasis):
	#spline basis(:,j) vs z at z
	if INTERP_BASIS_ORDER == 1:
		basis_j = interpolate.interp1d(zs, basis[:,j+1], kind='linear')
	elif INTERP_BASIS_ORDER == 3:
		basis_j = interpolate.interp1d(zs, basis[:,j+1], kind='cubic')
	basis_list.append(basis_j)	
		
def basis_eval(z, reion_nbasis):
	if ((z > zmax_basis) or (z < zmin)):	
		basis_val = np.zeros(reion_nbasis)
	else:
		basis_val = np.empty(reion_nbasis)	
		for j in np.arange(0, reion_nbasis):
			basis_val[j] = basis_list[j](z)
	return basis_val
	
#----------------------------------------------------------------  

helium_fullreion_redshift = 3.5
helium_fullreion_deltaredshift = 0.5
helium_fullreion_redshiftstart = 5.0
fHe = yheused/(mass_ratio_He_H * (1.0 - yheused))
print 'fHe = ', fHe, '; expecting ', 0.818405673350E-01

def xe_helium_eval(z):
	xod = (1.0+helium_fullreion_redshift - (1.0+z))/helium_fullreion_deltaredshift
	if (xod > 100.0):
		th = 1.0
	else:
		th = np.tanh(xod)
	xe_helium_val = fHe * (1.0 + th)/2.0	
	return xe_helium_val	

if LOAD_XEZ_RECOMB:
	# 5pc xe of recomb from camb for testing: z = 3 to 40
	xez_recomb = np.genfromtxt(fn_xez_recomb) # 2 cols: z vs xe(z) (used in this code)
	print 'recomb xez has z range: %f to %f'%(xez_recomb[0,0],xez_recomb[-1,0])			
	xe_recomb_interpolate = interpolate.interp1d(xez_recomb[:,0], xez_recomb[:,1], kind='linear')
	xe_recomb_eval = np.vectorize(xe_recomb_interpolate)

def xe_fid_eval(z): # include xe_fid and helium for now; might include recom and Gauss smoothing
	xlowz = 1.0 + fHe # or fraction
	#xe_recomb = 0.298082922827E-03 # 5pc, might need to generalize later
	#xe_recomb = 0.0
	if z <= zmin_basis: # zmin_basis
		xe_fid_val = xlowz 
	elif z >= zmax_basis: # zmax_basis
		#xe_fid_val = xe_recomb # one number for now, do interpolation later
		xe_fid_val = xe_recomb_eval(z) # one number for now, do interpolation later
	else:
		if z < zs[1]: 
			xe_fid_val = xlowz + (xe_fid - xlowz)/dz_basis * (z - zmin_basis)
			# valid for helium_fullreion_redshiftstart < zmin_basis
		elif z > zs[-2]:
			xe_recomb_tmp = xe_recomb_eval(zmax_basis)
			xe_fid_val = xe_fid + (xe_recomb_tmp - xe_fid)/dz_basis * (z - zs[-2])
		else:
			xe_fid_val = xe_fid		
	if z < helium_fullreion_redshiftstart:
		xe_fid_val = xe_fid_val + xe_helium_eval(z)
		#if z < zs[-2]:	
		#	xe_fid_val  = xe_fid_val - xe_recomb_eval(z)	
	return xe_fid_val 

class XeSmooth():

	def __init__(self, dz_fine_center, dz_fine_range_left, dz_fine_range_right, dz_fine, sigma_smoothing, reion_nbasis):
		self.dz_fine_center = dz_fine_center
		self.dz_fine_range_left = dz_fine_range_left
		self.dz_fine_range_right = dz_fine_range_right
		self.dz_fine = dz_fine
		self.sigma_smoothing= sigma_smoothing
		self.reion_nbasis = reion_nbasis
		[self.fine_grid, self.n_grid_points] = self.make_fine_grid()
		self.basis_fine = self.make_basis_fine()
		
	def make_fine_grid(self):
		fine_grid = np.arange(self.dz_fine_center - self.dz_fine_range_left, self.dz_fine_center + self.dz_fine_range_right + self.dz_fine, self.dz_fine)
		n_grid_points = fine_grid.size
		return fine_grid, n_grid_points
		
	def make_basis_fine(self):
		basis_fine = np.empty((self.n_grid_points, self.reion_nbasis))	
		for i in np.arange(0, self.n_grid_points):
			basis_fine[i,:] = basis_eval(self.fine_grid[i], reion_nbasis) 
		return basis_fine	

	def __call__(self, z, reion_mj):
		sigma2 = 2.0*self.sigma_smoothing**2.0
		eval = 0.0
		norm = 0.0
		xe_temp = xe_fid_eval(self.fine_grid[0]) + np.dot(reion_mj, self.basis_fine[0])
		gauss = np.exp(-(np.log((1.0+self.fine_grid[0])/(1.0+z)))**2/sigma2)/(1.0+self.fine_grid[0])
		eval = eval + xe_temp * gauss
		norm = norm + gauss
		for i in np.arange(1, self.fine_grid.size-1):
			xe_temp = xe_fid_eval(self.fine_grid[i]) + np.dot(reion_mj, self.basis_fine[i])
			gauss = np.exp(-(np.log((1.0+self.fine_grid[i])/(1.0+z)))**2/sigma2)/(1.0+self.fine_grid[i])
			eval = eval + 2.0 * xe_temp * gauss
			norm = norm + 2.0 * gauss
		xe_temp = xe_fid_eval(self.fine_grid[-1]) + np.dot(reion_mj, self.basis_fine[-1])
		gauss = np.exp(-(np.log((1.0+self.fine_grid[-1])/(1.0+z)))**2/sigma2)/(1.0+self.fine_grid[-1])
		eval = eval + xe_temp * gauss
		norm = norm + gauss
		# final result
		eval = 0.5 * eval * self.dz_fine
		norm = 0.5 * norm * self.dz_fine
		eval = eval/norm
		return eval			
	
# zstart = (1+xe_basis%zmax)*exp(8.0*xe_basis%fine%sigma)-1
# if(log((1._dl+Reionization_maxz)/(1._dl+xe_basis%zmax)) & 
     #   <16.0_dl*xe_basis%fine%sigma) then
# reionization_maxz = (1+xe_basis%zmax)*exp(16.0*xe_basis%fine%sigma)-1
		
dz_apply_smoothing = 5.0	
zmax_basis = 30.0
zmin_basis = 6.0		
reion_nbasis = 5
def xe_pc_eval(z, reion_nbasis):
	if ((z >= zmax_basis - dz_apply_smoothing) & (z <= zmax_basis)):
		dz_fine_center = zmax_basis			
		dz_fine_range_left = 5.0
		dz_fine_range_right = 5.0
		dz_fine = dz_basis/7.0
		sigma_smoothing = 0.015
		xe_smoothed_highz_eval = XeSmooth(dz_fine_center, dz_fine_range_left, dz_fine_range_right, dz_fine, sigma_smoothing, reion_nbasis)
		xe = xe_smoothed_highz_eval(z, reion_mj)
	elif ((z <= zmin_basis + dz_apply_smoothing) & (z >= zmin_basis)):
		dz_fine_center = zmin_basis			
		dz_fine_range_left = 0.0
		dz_fine_range_right = 5.0
		dz_fine = dz_basis/7.0
		sigma_smoothing = 0.015
		xe_smoothed_lowz_eval = XeSmooth(dz_fine_center, dz_fine_range_left, dz_fine_range_right, dz_fine, sigma_smoothing, reion_nbasis)
		xe = xe_smoothed_lowz_eval(z, reion_mj)
	else:
		xe = xe_fid_eval(z) + np.dot(reion_mj, basis_eval(z, reion_nbasis))			
	if USE_XEZ_RECOMB:
		xe  = xe - xe_recomb_eval(z)	
	return xe

#reion_redshift = ...
#reion_delta_redshift = 0.5
reion_zexp = 1.5
WindowVarMid = 0.40684689150538E+02
WindowVarDelta = 0.25795162261908E+01
#WindowVarDelta = reion_zexp*(1.0+reion_redshift)**(reion_zexp-1.0)*reion_delta_redshift
#WindowVarMid = (1.0+reion_redshift)**reion_zexp
xstart = 0.20525912929507E-03
def xe_hydrogen_eval(z):
	xod = (WindowVarMid - (1.0+z)**reion_zexp)/WindowVarDelta
	if (xod > 100.0):
		tgh=1.0
	else:
		tgh=np.tanh(xod)
	return ((1.0+fHe) - xstart)*(tgh+1.0)/2.0 + xstart
	
def xe_tanh_eval(z):
	xe = xe_hydrogen_eval(z)
	if z < helium_fullreion_redshiftstart:
		xe = xe + xe_helium_eval(z)
	if USE_XEZ_RECOMB:
		xe = xe - xe_recomb_eval(z)	
	return xe	
		
#----------------------------------------------------------------  

# calculate densities and H
# grho = 8piG * rho/c^2
#grhom = 3.0*H0**2.0/c**2.0*1000.0**2.0 #3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
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

Nnow = omegab*(1.0-yheused)*grhom*c**2.0/kappa/m_H/Mpc**2.0
print 'Nnow = ', Nnow, '; expecting', 0.189222966061E+00 # 0.189225846211

akthom = sigma_thomson*Nnow*Mpc

# captured: grhob, grhoc, grhog, grhonomass, grhok, grhol
def dtda(z):
	a = 1.0/(1.0+z)
	a2 = a**2
	#  8*pi*G*rho*a**4.
	grhoa2 = grhok*a2 + (grhoc + grhob) *a + grhog + grhornomass + grhol*a2**2
	if (num_nu_massive != 0):
		grhoa2=grhoa2 + grhonu_mass_eval(z)
	dtda = np.sqrt(3.0/grhoa2)
	return dtda
dtda_vec = np.vectorize(dtda)

# dtauda from camb for testing: log in z from 4 to -13
#fn_dtauda = 'test8_5pc_dtauda.dat'
fn_dtauda = 'test8_tanh_base2p6_omega_vary1_dtauda.dat'
dtauda = np.genfromtxt(fn_dtauda) # 2 cols: z vs xe(z) (used in this code)
print 'dtauda has z range: %f to %f'%(dtauda[0,0],dtauda[-1,0])			
dtauda_interpolate = interpolate.interp1d(dtauda[:,0], dtauda[:,1], kind='cubic')
dtauda_eval = np.vectorize(dtauda_interpolate)


#if USE_XEZ_CMB:
# tanh xe(z) from camb for testing: z = 3 to 40
xez_camb = np.genfromtxt(fn_xez_camb) # 2 cols: z vs xe(z) (used in this code)
print 'camb xez has z range: %f to %f'%(xez_camb[0,0],xez_camb[-1,0])			
xe_camb_interpolate = interpolate.interp1d(xez_camb[:,0], xez_camb[:,1], kind='cubic')
xe_camb_eval = np.vectorize(xe_camb_interpolate)

# omeganuh2 evolution in z output from camb: z = 3 to 40
nu_evol = np.genfromtxt(fn_nu_evol) # 2 cols: z vs grhonu_mass(z) (used in this code)
print 'nu_evol has z range: %f to %f'%(nu_evol[0,0], nu_evol[-1,0])		
#omeganu_evol = nu_evol[:,1]*(1.0+nu_evol[:,0])**4/3.0/(H0**2.0)
#omeganu_evol_interpolate = interpolate.interp1d(nu_evol[:,0], omega_nu_evol, kind='cubic')
grhonu_mass_interpolate = interpolate.interp1d(nu_evol[:,0], nu_evol[:,1], kind='cubic')
grhonu_mass_eval = np.vectorize(grhonu_mass_interpolate)

#def square_ratio_H_H0_eval(z):
#	x = 1.0 + z
#	return omegag * x**4.0 + (omegab + omegam) * x**3.0 + omegal + nu_evol_eval(z) #k = 0 

#def xe_eval(z): make later with mj

def integrand(z):	
	#return xe_eval(z) * (1.0 + z)**2.0 * H0_over_H_eval(z)
	#return xe_camb_eval(z) * (1.0 + z)**2.0 * np.sqrt(square_ratio_H_H0_eval(z))
	if USE_XEZ_CMB:
		integ = -xe_camb_eval(z) * dtda(z)
		#integ = -xe_camb_eval(z) * dtauda_eval(z)
	else:
		integ = xe_pc_eval(z, reion_nbasis) * dtda(z) # pcs
		#integ = xe_tanh_eval(z) * dtda(z) # tanh
	return integ
		
def calculate_tau(z1, z2):
	print 'integration for tau from z1 = ', z1, ' to z2 = ', z2
	limit_nstep = 1000
	int = integrate.quad(integrand, z1, z2, limit=limit_nstep)
	print int
	tau = akthom * int[0]
	return tau

def calculate_tau_total(tau, tau_zmin):
	return tau + tau_zmin
	
# Report result of tau(z1, z2) and tau_total
tau = calculate_tau(z1, z2)
print 'tau = ', tau
print 'tau_total  = ', tau + tau_zmin + tau_zmax
#if z1 == zmin:
#	tau_total = calculate_tau_total(tau, tau_zmin)
#	print 'tau_total  = ', tau_total

#dtauda = np.genfromtxt(fn_dtauda) 
#dtauda_me = dtda_vec(dtauda[5:800,0])

#print (dtauda_me - dtauda)/dtauda


fig, ax = plt.subplots(1,figsize=(6,4.5))
xe_pc = np.empty(zs.size)
i = 0
for z in zs:
	xe_pc[i] = xe_pc_eval(z, reion_nbasis)
	i = i + 1 
plt.plot(zs, xe_pc, 'r-')
plt.plot(zs, -xe_camb_eval(zs), 'b-')

legloc1 = 'upper right'
title1 = r'$x_e(z)$ python vs cosmomc for 5pcs'
xlabel1 = r'$z$'
ylabel1 = r'$x_e$'

plt.legend(loc=legloc1)	
plt.title(title1)
ax = plt.gca()
ax.set_xlabel(xlabel1)	
ax.set_ylabel(ylabel1)	
ax.set_xlim(np.array([z1, z2]))

