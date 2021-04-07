import numpy as np
import scipy.interpolate as interpolate
import ReionBasis

# ------------------------------- Settings ----------------------------

# files: xe_recomb(z), basis
# internal keys: is_mortonson_format, INTERP_KIND_XEZ_CAMB
# might want to make global keys: USE_XEZ_RECOMB, xe_fid, yheused

fn_xez_recomb = 'python/tanh_base2p6_xez_recomb.dat' # cosmomc xe(z) recomb used for testing
fn_basis = 'python/xepcs.dat' # not original mortonson file
is_mortonson_format = True # only support true right now

Basis = ReionBasis.ReionBasis(fn_basis, is_mortonson_format)

USE_XEZ_RECOMB = True
INTERP_KIND_XEZ_CAMB = 'cubic'

xe_fid = 0.15
helium_fullreion_redshift = 3.5
helium_fullreion_deltaredshift = 0.5
helium_fullreion_redshiftstart = 5.0
yheused = 0.245341   # TO DO: need to pass this in 
mass_ratio_He_H =   0.397150000000E+01
fHe = yheused/(mass_ratio_He_H * (1.0 - yheused))
print 'fHe = ', fHe, '; expecting ', 0.818405673350E-01


# ------------------------ functions and class XeSmooth ---------------------------

def xe_helium_eval(z):
	xod = (1.0+helium_fullreion_redshift - (1.0+z))/helium_fullreion_deltaredshift
	if (xod > 100.0):
		th = 1.0
	else:
		th = np.tanh(xod)
	xe_helium_val = fHe * (1.0 + th)/2.0	
	return xe_helium_val	


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


def xe_fid_eval(z): # include xe_fid and helium for now; might include recom and Gauss smoothing
	xlowz = 1.0 + fHe # or fraction
	#xe_recomb = 0.298082922827E-03 # 5pc, might need to generalize later
	xe_recomb = 0.0
	if z <= Basis.zs[0]: # zmin_basis
		xe_fid_val = xlowz 
	#elif z >= Basis.zs[-1]: # zmax_basis
		#xe_fid_val = xe_recomb # one number for now, do interpolation later
	else:
		if z < Basis.zs[1]: 
			xe_fid_val = xlowz + (xe_fid - xlowz)/Basis.dz_basis * (z - Basis.zmin_basis)
			# valid for helium_fullreion_redshiftstart < zmin_basis
		#elif z > Basis.zs[-2]:
		#	xe_fid_val = xe_fid + (xe_recomb - xe_fid)/Basis.dz_basis * (z - Basis.zmax_basis)
		else:
			xe_fid_val = xe_fid		
	if z < helium_fullreion_redshiftstart:
		xe_fid_val = xe_fid_val + xe_helium_eval(z)
		#if z < Basis.zs[-2]:	
		#	xe_fid_val  = xe_fid_val - xe_recomb_eval(z)	
	return xe_fid_val 

def xe_pc_eval(z, reion_nbasis, reion_mj):
	xe_fid_fcn = xe_fid_eval(z)
	xe = xe_fid_fcn + np.dot(reion_mj, Basis.BasisEval(reion_nbasis, z))
	if USE_XEZ_RECOMB:
		xe  = xe - xe_recomb_eval(z)	
	return xe
	
if USE_XEZ_RECOMB:
	# 5pc xe of recomb from camb for testing: z = 3 to 40
	xez_recomb = np.genfromtxt(fn_xez_recomb) # 2 cols: z vs xe(z) (used in this code)
	print 'recomb xez has z range: %f to %f'%(xez_recomb[0,0],xez_recomb[-1,0])			
	xe_recomb_interpolate = interpolate.interp1d(xez_recomb[:,0], xez_recomb[:,1], kind=INTERP_KIND_XEZ_CAMB)
	xe_recomb_eval = np.vectorize(xe_recomb_interpolate)


#-------------------------- For tanh reionization ----------------------------  

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
	
# get tanh xe(z) from camb for testing: z = 3 to 40
#fn_xez_camb = 'test8_tanh_xez_with_helium_without_recomb.dat' # cosmomc xe(z) used for testing
#xez_camb = np.genfromtxt(fn_xez_camb) # 2 cols: z vs xe(z) (used in this code)
#print 'camb xez has z range: %f to %f'%(xez_camb[0,0],xez_camb[-1,0])			
#xe_camb_interpolate = interpolate.interp1d(xez_camb[:,0], xez_camb[:,1], kind='cubic')
#xe_camb_eval = np.vectorize(xe_camb_interpolate)
