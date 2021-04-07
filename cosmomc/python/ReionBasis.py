import numpy as np
import scipy.interpolate as interpolate

class ReionBasis:

	# varaible shared by all instances
	
	def __init__(self, fn_basis, is_mortonson_format):
		self.reion_nbasis_max = 10 # boost this number if call BasisEval with reion_nbasis > self.reion_nbasis_max
		self.fn_basis = fn_basis
		self.is_mortonson_format = is_mortonson_format
		self.basis_read_done = False
		self.basis_setup_done = False
		self.BasisInit()
			
	def BasisRead(self): # only support is_mortonson_format = True so far, need to test False case.
		fn_basis = self.fn_basis
		is_mortonson_format = self.is_mortonson_format
		if is_mortonson_format:
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
		dz_basis = zs[1] - zs[0] # assuming increasing z order ! CAN CHECK
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
		
		basis_read_done = True
		
		return zs, basis, nz, nbasis, zmin_basis, zmax_basis, dz_basis, basis_read_done
		
	def BasisSetup(self):
		basis_list = []
		for j in np.arange(0, self.reion_nbasis_max):
			#spline basis(:,j) vs z at z
			basis_j = interpolate.interp1d(self.zs, self.basis[:,j+1], kind='linear')
			basis_list.append(basis_j)	
			basis_setup_done = True
		return basis_list, basis_setup_done	
				
	def BasisInit(self):
		[self.zs, self.basis, self.nz, self.nbasis, self.zmin_basis, self.zmax_basis, self.dz_basis, self.basis_read_done] = self.BasisRead()
		[self.basis_list, self.basis_setup_done] = self.BasisSetup()
		print 'Done reading and setting up basis.'
	
	def BasisEval(self, reion_nbasis, z):
		if reion_nbasis > self.reion_nbasis_max:
			print 'Error from BasisEval: you asked for reion_nbasis = %i > reion_nbasis_max = %i. Increase reion_nbasis_max in ReionBasis.py.'%(reion_nbasis, self.reion_nbasis_max)
		else:
			basis_val = np.empty(reion_nbasis)
			for j in np.arange(0, reion_nbasis):
				basis_val[j] = self.basis_list[j](z)
			return basis_val
				