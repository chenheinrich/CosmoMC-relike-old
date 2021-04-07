import getdist.plots as plots, os
g=plots.GetDistPlotter(plot_data=r'./plot_data/test6/')
g.settings.setWithSubplotSize(4.0)
roots = ['test6']
g.triangle_plot(roots, ['reion_m1', 'reion_m2', 'reion_m3', 'reion_m4', 'reion_m5'])
g.export(os.path.join(r'results/test6/',r'test6_tri.pdf'))
