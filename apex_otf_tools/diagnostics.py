import pylab as pl
import mpl_plot_templates
import numpy as np

def diagplot(data, tsys, noise, dataset):

    pl.figure(1)
    pl.clf()
    mpl_plot_templates.imdiagnostics(data)
    pl.savefig(dataset+"_diagnostics.pdf",bbox_inches='tight')
    pl.figure(2)
    pl.clf()
    pl.subplot(2,1,1)
    pl.plot(tsys,np.arange(tsys.size),alpha=0.5)
    pl.xlabel("TSYS")
    pl.ylabel("Integration")
    pl.subplot(2,1,2)
    pl.plot(tsys, noise, '.',alpha=0.5)
    pl.xlabel("TSYS")
    pl.ylabel("Noise")
    pl.savefig(dataset+"_tsys.pdf",bbox_inches='tight')
