import getdist.plots as plots, os
g=plots.GetDistPlotter(plot_data=r'./plot_data/')
g.settings.setWithSubplotSize(4.0)
roots = ['test2']
markers={}
g.plots_1d(roots, markers=markers)
g.export(os.path.join(r'',r'test2.pdf'))
