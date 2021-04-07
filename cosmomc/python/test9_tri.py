import getdist.plots as plots, os
g=plots.GetDistPlotter(plot_data=r'./plot_data/test9/')
g.settings.setWithSubplotSize(4.0)
roots = ['test9']
g.triangle_plot(roots, ['omegabh2', 'omegach2', 'tau', 'ns', 'logA', 'H0', 'omegam', 'omegal', 'sigma8'])
g.export(os.path.join(r'results/test9/',r'test9_tri.pdf'))
