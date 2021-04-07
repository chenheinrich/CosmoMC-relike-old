import getdist.plots as gplot
g = gplot.getSinglePlotter(chain_dir=r'chains/cosmomc_test2_chains/')
#g = gplot.getSinglePlotter(chain_dir=r'chains/')

samples = g.sampleAnalyser.samplesForRoot('test2')
#samples = g.sampleAnalyser.samplesForRoot('test6')

p = samples.getParams()

derived = p.sigma8 * p.omegam ** 0.6

print 'mean, err = ', samples.mean(derived), samples.std(derived)
print '95% limits: ', samples.twoTailLimits(derived, 0.95)

#summj = p.reion_m1 + p.reion_m2 + p.reion_m3 + p.reion_m4 + p.reion_m5
#print 'mean, err = ', samples.mean(summj), samples.std(summj)
#print '95% limits: ', samples.twoTailLimits(summj, 0.95)

rd_fid = 149.28
rsH = p.Hubble057 * p.rdrag / rd_fid
samples.addDerived(rsH, name='rsH', label=r'H(0.57) (r_{\mathrm{drag}}/r_{\mathrm{drag}}^{\rm fid})\, [{\rm km} \,{\rm s}^{-1}{\rm Mpc}^{-1}]')
samples.updateBaseStatistics()
