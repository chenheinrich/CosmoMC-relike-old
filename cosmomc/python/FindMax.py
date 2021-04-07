import numpy as np
from matplotlib import pyplot as plt

chain_dir = '/project/kicp/chenhe/Reionization/cosmomc/chains/test6_chains/'
root = 'test6'
i = 1
sample = np.genfromtxt(chain_dir + root + '_%s'%i + '.txt')
for i in np.arange(2,5):
	sample_tmp = np.genfromtxt(chain_dir + root + '_%s'%i + '.txt')
	sample = np.vstack((sample, sample_tmp))
	
# mjs: cols (5 to 9) + 2, index (4 to 8) + 2
sample_like = sample[:,1] 
indp_m1 = 4
reion_nbasis = 5
indc_m1 = indp_m1+2
m1 = sample[:,indc_m1 ] 
m2 = sample[:,indc_m1 +1] 	
m3 = sample[:,indc_m1 +2] 
m4 = sample[:,indc_m1 +3] 
m5 = sample[:,indc_m1 +4] 
#print sample_like[0:5]

ms_sample = sample[:,indc_m1:indc_m1+reion_nbasis]

# find and print results for large m5 > 0.4, max likelihood point

ind_m5 = np.where(m5>0.4)
like_m5 = sample_like[ind_m5]
max_ind_like = like_m5.argmin()
max_like = like_m5[max_ind_like]
max_ind_m5_like = ind_m5[0][max_ind_like] # ind_m5 is a tuple # 15332

print 'model to test (large m5 > 0.4, max likelihood point)'
print 'index = ', max_ind_m5_like  # 126448
print 'm1 = ', m1[max_ind_m5_like] # 0.06909047
print 'm2 = ', m2[max_ind_m5_like] # 0.1085179
print 'm3 = ', m3[max_ind_m5_like] # 0.09189145
print 'm4 = ', m4[max_ind_m5_like] # -0.1108808
print 'm5 = ', m5[max_ind_m5_like] # 0.4557212
print 'large m5 > 0.4, max likelihood point'
print sample[max_ind_m5_like, :]
print 'likelihood = ', sample_like[max_ind_m5_like] # 
print 'global max likelihood = ', np.min(sample_like) # 
like_lowTEB = sample[max_ind_m5_like,-4]/2 
like_plik = sample[max_ind_m5_like,-3]/2 
like_prior = sample[max_ind_m5_like,-2]/2 
like_CMB = sample[max_ind_m5_like,-1]/2 

print 'like_lowTEB = ', like_lowTEB # 5248.92
print 'like_plik = ', like_plik # 1239.3745
print 'like_prior = ', like_prior # 10.74479
print 'like_CMB = ', like_CMB # 6488.295

# find and print results for global max likelihood point

gmax_ind_like = sample_like.argmin() 
print 'model to test (global max likelihood point)'
print 'index = ', gmax_ind_like # 100578 for test 6
print 'm1 = ', m1[gmax_ind_like]
print 'm2 = ', m2[gmax_ind_like]
print 'm3 = ', m3[gmax_ind_like]
print 'm4 = ', m4[gmax_ind_like]
print 'm5 = ', m5[gmax_ind_like]
print 'global max likelihood point:'
print sample[gmax_ind_like, :]
print 'likelihood = ', sample_like[gmax_ind_like]

like_lowTEB = sample[gmax_ind_like,-4]/2
like_plik = sample[gmax_ind_like,-3]/2
like_prior = sample[gmax_ind_like,-2]/2
like_CMB = sample[gmax_ind_like,-1]/2

print 'like_lowTEB = ', like_lowTEB
print 'like_plik = ', like_plik
print 'like_prior = ', like_prior
print 'like_CMB = ', like_CMB

# ----- Check parameter difference (Cls differ by 4% at high l)

As1 =  np.min(sample[gmax_ind_like,11])
As2 = np.min(sample[max_ind_m5_like,11])
As1_true = np.exp(As1)
As2_true = np.exp(As2)
(As1_true - As2_true)/As2_true # -0.040732205886480459
# As differ by 4% as well.


# ----- Print foreground parameters for maxlike (for standalone like) ------

sample_fg = sample[gmax_ind_like, 13:40]
acib217 = sample_fg[1]
ncib = -1.3
galfEEindex = -2.4
galfTEindex = -2.4
calpol = 1.0
calEE0 = 1.0
calEE1 = 1.0
calEE2 = 1.0
A_planck = sample_fg[0]
cal0 = sample_fg[25]
cal2 = sample_fg[26]

fg = np.zeros(94)
fg[0] = acib217 #A_cib
fg[1] = ncib # ncib or cib_index
fg[2:19] = sample_fg[2:19]
fg[19] = galfEEindex # galf_EE_index
fg[20:26] = sample_fg[19:25]
fg[26] = galfTEindex # galf_TE_index
fg[-7] = cal0 # calib_100T
fg[-6] = cal2 # calib_217T
fg[-5] = calEE0 # calib_100P
fg[-4] = calEE1 # calib_143P
fg[-3] = calEE2 # calib_217P
fg[-2] = calpol # A_pol
fg[-1] = A_planck

fn_fg = '/home/chenhe/Projects/Reionization/reionization2016/cosmomc/test6_maxlike_hi_full.dat'
np.savetxt(fn_fg, fg)
print 'file saved: ', fn_fg


# ----- Print foreground parameters for m5>0.4 max likelihood(for standalone like) ------

sample_fg = sample[max_ind_m5_like, 13:40]
acib217 = sample_fg[1]
ncib = -1.3
galfEEindex = -2.4
galfTEindex = -2.4
calpol = 1.0
calEE0 = 1.0
calEE1 = 1.0
calEE2 = 1.0
A_planck = sample_fg[0]
cal0 = sample_fg[25]
cal2 = sample_fg[26]

fg = np.zeros(94)
fg[0] = acib217 #A_cib
fg[1] = ncib # ncib or cib_index
fg[2:19] = sample_fg[2:19]
fg[19] = galfEEindex # galf_EE_index
fg[20:26] = sample_fg[19:25]
fg[26] = galfTEindex # galf_TE_index
fg[-7] = cal0 # calib_100T
fg[-6] = cal2 # calib_217T
fg[-5] = calEE0 # calib_100P
fg[-4] = calEE1 # calib_143P
fg[-3] = calEE2 # calib_217P
fg[-2] = calpol # A_pol
fg[-1] = A_planck

fn_fg = '/home/chenhe/Projects/Reionization/reionization2016/cosmomc/test6_largem5_hi_full.dat'
np.savetxt(fn_fg, fg)
print 'file saved: ', fn_fg



#0 calPlanck	y_{\rm cal}
#1 acib217	A^{CIB}_{217}

#2 xi	\xi^{tSZ-CIB}
#3 asz143	A^{tSZ}_{143}
#4 aps100	A^{PS}_{100}
#5 aps143	A^{PS}_{143}
#6 aps143217	A^{PS}_{143\times217}
#7 aps217	A^{PS}_{217}
#8 aksz	A^{kSZ}
#9 kgal100	A^{{\rm dust}TT}_{100}
#10 kgal143	A^{{\rm dust}TT}_{143}
#11 kgal143217	A^{{\rm dust}TT}_{143\times217}
#12 kgal217	A^{{\rm dust}TT}_{217}
#13 galfEE100	A^{{\rm dust}EE}_{100}
#14 galfEE100143	A^{{\rm dust}EE}_{100\times143}
#15 galfEE100217	A^{{\rm dust}EE}_{100\times217}
#16 galfEE143	A^{{\rm dust}EE}_{143}
#17 galfEE143217	A^{{\rm dust}EE}_{143\times217}
#18 galfEE217	A^{{\rm dust}EE}_{217}

#19 galfTE100	A^{{\rm dust}TE}_{100}
#20 galfTE100143	A^{{\rm dust}TE}_{100\times143}
#21 galfTE100217	A^{{\rm dust}TE}_{100\times217}
#22 galfTE143	A^{{\rm dust}TE}_{143}
#23 galfTE143217	A^{{\rm dust}TE}_{143\times217}
#24 galfTE217	A^{{\rm dust}TE}_{217}

#25 cal0	c_{100}
#26 cal2	c_{217}
