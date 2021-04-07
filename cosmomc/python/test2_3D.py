import getdist.plots as plots, os
g=plots.GetDistPlotter(plot_data=r'./plot_data/')
g.settings.setWithSubplotSize(6.0)
roots = ['test2']
sets=[]
sets.append(['H0','omegam','tau'])
g.plots_3d(roots,sets)
g.export(os.path.join(r'',r'test2_3D.pdf'))
