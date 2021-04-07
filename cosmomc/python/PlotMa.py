import numpy as np
from matplotlib import pyplot as plt

import getdist.plots as gplot
g = gplot.getSinglePlotter(chain_dir=r'/project/kicp/chenhe/Reionization/cosmomc/chains/')
  
samples = g.sampleAnalyser.samplesForRoot('test6')
p = samples.getParams()
p.reion_m1

pair_2D = np.array([[3,5]])
for i, j in pair_2D:
	roots = ['test6']
	g = gplot.getSinglePlotter(chain_dir=r'/project/kicp/chenhe/Reionization/cosmomc/chains/')
	g.plot_2d(roots, 'reion_m%s'%i, 'reion_m%s'%j, filled=True)
	g.add_legend(['plik full'], legend_loc='upper right');
	g.export('test6_2D_m%s_m%s.pdf'%(i,j))

g = gplot.getSubplotPlotter(chain_dir=r'/project/kicp/chenhe/Reionization/cosmomc/chains/')
roots = ['test6']
reion_nbasis = 5
mvars = ["reion_m%s"%i for i in np.arange(1,reion_nbasis+1)]
g.triangle_plot(roots, mvars, filled=True)
g.export('test6_mj_triangle.pdf')

for i in np.arange(1,reion_nbasis+1):
	g = gplot.getSinglePlotter(chain_dir=r'/project/kicp/chenhe/Reionization/cosmomc/chains/', width_inch=4)
	g.plot_1d(samples, 'reion_m%s'%i)
	g.export('test6_1D_m%s.pdf'%i)
	
from getdist.densities import Density2D
import getdist.plots as gplot
import numpy as np

g = gplot.getSinglePlotter(chain_dir=r'./PLA')
xvalues = np.arange(85, 110, 0.3)
yvalues = np.arange(1000, 1500, 4)
x,y = np.meshgrid(xvalues, yvalues)
loglike = my_loglike_func(x,y)
density = Density2D(xvalues,yvalues, np.exp(-loglike / 2))
density.contours = np.exp(-np.array([1.509, 2.4477]) ** 2 / 2)
g.add_2d_contours(root, 'x', 'y', filled=True, density=density)
	