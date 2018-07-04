from __future__ import print_function, division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
rc('font', family='serif')
from operator import itemgetter
import os
import math

from astropy.table import Table
from astropy.io import fits, ascii
from astropy import units as u

AP_PATH = os.environ['APOGEE_DATA']

import apogee_tools as ap


def plotCrossValidation(trn_labels, crv_labels, **kwargs):

    # required
    label_names = kwargs.get('label_names', ['TEFF', '[Fe/H]'])

    # optional
    save = kwargs.get('save', False)
    out  = kwargs.get('out', 'Cross_Validation_Scatter.pdf')

    if len(label_names) == 1:
        fig, ax1 = plt.subplots(1, 1)

        trn_teffs = list(map(itemgetter(0), trn_labels))

        crv_teffs = list(map(itemgetter(0), crv_labels))

        ax1.scatter(trn_teffs, crv_teffs, color='k')

        limits_1 = [ np.min([ax1.get_xlim(), ax1.get_ylim()]), \
                 np.max([ax1.get_xlim(), ax1.get_ylim()]) ]

        # now plot both limits against eachother
        ax1.plot(limits_1, limits_1, '--', alpha=0.75, color='r', zorder=0, linewidth=2)
        ax1.set_aspect('equal')
        ax1.set_xlim(limits_1)
        ax1.set_ylim(limits_1)

        ax1.set_title(label_names[0] + 'Labels')
        ax1.set_xlabel('Training Label')
        ax1.set_ylabel('Cross-Validated Label')

    elif len(label_names) == 2:
        fig, (ax1, ax2) = plt.subplots(1, 2)

        trn_teffs = list(map(itemgetter(0), trn_labels))
        trn_metal = list(map(itemgetter(1), trn_labels))

        crv_teffs = list(map(itemgetter(0), crv_labels))
        crv_metal = list(map(itemgetter(1), crv_labels))

        ax1.scatter(trn_teffs, crv_teffs, color='k')

        limits_1 = [ np.min([ax1.get_xlim(), ax1.get_ylim()]), \
                 np.max([ax1.get_xlim(), ax1.get_ylim()]) ]

        # now plot both limits against eachother
        ax1.plot(limits_1, limits_1, '--', alpha=0.75, color='r', zorder=0, linewidth=2)
        ax1.set_aspect('equal')
        ax1.set_xlim(limits_1)
        ax1.set_ylim(limits_1)

        ax1.set_title(label_names[0] + ' Labels')
        ax1.set_xlabel('Training Label')
        ax1.set_ylabel('Cross-Validated Label')

        ax2.scatter(trn_metal, crv_metal, color='k')

        limits_2 = [ np.min([ax2.get_xlim(), ax2.get_ylim()]), \
                 np.max([ax2.get_xlim(), ax2.get_ylim()]) ]

        ax2.plot(limits_2, limits_2, '--', alpha=0.75, color='r', zorder=0, linewidth=2)
        ax2.set_aspect('equal')
        ax2.set_xlim(limits_2)
        ax2.set_ylim(limits_2)

        ax2.set_title(label_names[1] + ' Labels')
        ax2.set_xlabel('Training Label')
        # ax2.set_ylabel('Cross-Validated Label')

    plt.savefig(str(out))

    plt.tight_layout()
    plt.show()
    plt.close()


def plotCannonModels(ds, te_flux, te_labels, **kwargs):

    bands  = [[15160,15800],[15880,16420],[16500,16935]]
    n      = kwargs.get('band', 1)
    yrange = kwargs.get('yrange', [.6,1.2])
    base   = kwargs.get('base', 'Ref')
    lbl_names = kwargs.get('lbl_names', ['SPT'])
    sigfig = kwargs.get('sigfig', [1 for x in lbl_names])
    snr = kwargs.get('snr', [])
    
    save = kwargs.get('save', False)
    out  = kwargs.get('out', 'Models_Band'+str(n)+'.pdf')
    
    nspecs = len(te_flux)
    nparam = len(lbl_names)

    nplots = kwargs.get('nplots', nspecs)
    
    tr_label = ds.tr_label
    tr_label_unc = kwargs.get('tr_lbl_unc')
    te_label = te_labels
    te_label_unc = kwargs.get('te_lbl_unc', [0,0])
    
    wl = ds.wl
    tr_flux = ds.tr_flux
    tr_ivar = ds.tr_ivar
    
    fig, axs = plt.subplots(nplots, 1, figsize=(12,3*nspecs))
    for i, ax in enumerate(fig.axes):

        tr_stdev = [1/math.sqrt(ivar) for ivar in tr_ivar[i]]
    
        data = ap.Spectrum(wave=wl, flux=tr_flux[i], sigmas=tr_stdev)
        mdl  = ap.Spectrum(wave=wl, flux=te_flux[i])
        chi  = ap.compareSpectra(data, mdl, fit_scale=False)[0]
        
        ax.plot(wl, tr_flux[i], color='k', linewidth=1)
        ax.plot(wl, te_flux[i], color='r', linewidth=1)
        
        par_str1 = [r'${} = {} \pm {}$'.format(lbl_names[k], round(te_label[i][k],sigfig[k]), te_label_unc[0]) for k in range(nparam)]
        can_lbl  = r'$Cannon: $ {}'.format(', '.join(par_str1))
        ax.text(bands[n-1][0]+10, yrange[1]-.08, can_lbl, color='r', fontsize=15, va='bottom', ha='left')

        par_str2 = [r'${} = {} \pm {}$'.format(lbl_names[k], round(tr_label[i][k],sigfig[k]), tr_label_unc[0]) for k in range(nparam)]
        ref_lbl  = r'${}: $ {}'.format(base, ', '.join(par_str2))
        ax.text(bands[n-1][0]+10, yrange[0]+.08, ref_lbl, color='k', fontsize=15, va='top', ha='left')
        
        chi_lbl = r'$\chi^{2} = %s$'%(str(chi))
        ax.text(bands[n-1][1]-10, yrange[1]-.08, chi_lbl, color='r', fontsize=15, va='bottom', ha='right')
        
        if len(snr) != 0:
            snr_lbl = r'$SNR = {}$'.format(str(round(snr[i],1)))
            ax.text(bands[n-1][1]-10, yrange[0]+.08, snr_lbl, color='k', fontsize=15, va='top', ha='right')
         
        ax.set_title(r'${}$'.format(ds.tr_ID[i]), fontsize=20)
        ax.set_xlim(bands[n-1])
        ax.set_ylim(yrange)
        ax.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]', fontsize=15)
        if i == nspecs-1:
            ax.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    
    plt.tight_layout()
    if save == True:
        plt.savefig(str(out))
    plt.show()
    plt.close()


def plotOHBands(**kwargs):

    """
    Plot OH-bands which are sensitive to Teff.
    Unless specified, 1) black = spectrum, 2) red = cannon model, 3) blue = btsettl, 4) green = other
    """

    # required
    specs = kwargs.get('specs')

    # optional 
    labels = kwargs.get('labels', ['Data', 'Cannon', 'BTSettl'])

    bands = [[15400,15450], [16350,16360], [16860,16890]]
    nbands = len(bands)

    fig, axs = plt.subplots(1, nbands, figsize=(12,3))

    nspecs = len(specs)
    colors = ['k', 'r', 'b', 'g']

    for i, ax in enumerate(fig.axes):
        for j in range(nspecs):
            ax.plot(specs[j].wave, specs[j].flux, color=colors[j])
        ax.set_xlim(bands[i])

    plt.tight_layout()
    plt.show()
    plt.close()


def plotSpectralSequence(**kwargs):

    """
    Plot spectral sequence, with flux offset. Optional: overplot models.
    """

    # required
    specs  = kwargs.get('specs')
    models = kwargs.get('models')
    labels = kwargs.get('labels')

    # optional
    save   = kwargs.get('save', True)
    out    = kwargs.get('out', 'Spectral_Sequence.pdf')

    nspecs = len(specs)

    plotbands = [[15200,15800],[15870,16420],[16500,16930]] 

    highlights1 = [[15290,15350],[15530,15590],[15620,15650],[15715,15780]] 
    highlights2 = [[15960,15980],[16150,16220]]
    highlights3 = [[16520,16560],[16695,16775]]
    
    fig1, ax1 = plt.subplots(1, 1, figsize=(16,8))

    for i in range(nspecs):
        # Replace zeroes with nans
        spec_fl = specs[i].flux
        spec_fl[spec_fl == 0] = np.nan
        mdl_fl = models[i].flux
        mdl_fl[mdl_fl == 0] = np.nan

        ax1.plot(specs[i].wave, spec_fl + .5*i, color='k', alpha=.8, label=labels[i])
        ax1.plot(models[i].wave, models[i].flux + .5*i, color='r', alpha=.8)
        ax1.text(plotbands[0][1]+10, .5*i + .9, 'M'+str(i), fontweight='bold', fontsize=16)
        ax1.text(plotbands[0][0]+10, .5*i + 1.1, str(specs[i].name), fontweight='bold', fontsize=12)
    for h in highlights1:
        ax1.axvspan(h[0], h[1], color='b', alpha=0.1)
    ax1.set_xlim(plotbands[0])
    ax1.set_ylim([.5,6])
    ax1.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    ax1.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$] + offset', fontsize=15)
    ax1.set_title('APOGEE M dwarf Spectral Sequence', fontsize=20)

    fig1.savefig('1_' + str(out))
    plt.show()
    plt.close(fig1)

    fig2, ax2 = plt.subplots(1, 1, figsize=(16,8))

    for i in range(nspecs):
        spec_fl = specs[i].flux
        spec_fl[spec_fl == 0] = np.nan
        mdl_fl = models[i].flux
        mdl_fl[mdl_fl == 0] = np.nan

        ax2.plot(specs[i].wave, spec_fl + .5*i, color='k', alpha=.8, label=labels[i])
        ax2.plot(models[i].wave, mdl_fl + .5*i, color='r', alpha=.8)
        ax2.text(plotbands[1][1]+10, .5*i + .9, 'M'+str(i), fontweight='bold', fontsize=16)
    for h in highlights2:
        ax2.axvspan(h[0], h[1], color='b', alpha=0.1)
    ax2.set_xlim(plotbands[1])
    ax2.set_ylim([.5,6])
    ax2.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    ax2.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$] + offset', fontsize=15)

    fig2.savefig('2_' + str(out))
    plt.show()
    plt.close(fig2)

    fig3, ax3 = plt.subplots(1, 1, figsize=(16,8))

    for i in range(nspecs):
        spec_fl = specs[i].flux
        spec_fl[spec_fl == 0] = np.nan
        mdl_fl = models[i].flux
        mdl_fl[mdl_fl == 0] = np.nan

        ax3.plot(specs[i].wave, spec_fl + .5*i, color='k', alpha=.8, label=labels[i])
        ax3.plot(models[i].wave, mdl_fl + .5*i, color='r', alpha=.8)
        ax3.text(plotbands[2][1]+10, .5*i + .9, 'M'+str(i), fontweight='bold', fontsize=16)
    for h in highlights3:
        ax3.axvspan(h[0], h[1], color='b', alpha=0.1)
    ax3.set_xlim(plotbands[2])
    ax3.set_ylim([.5,6])
    ax3.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=15)
    ax3.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$] + offset', fontsize=15)

    fig3.savefig('3_' + str(out))
    plt.show()
    plt.close(fig3)


def plotBands(**kwargs):

    """
    Plot OH-bands which are sensitive to Teff.
    Unless specified, 1) black = spectrum, 2) red = cannon model, 3) blue = btsettl, 4) green = other
    """

    # required
    specs = kwargs.get('specs')
    plot_band = kwargs.get('bands', 'OH')

    # optional 
    title  = kwargs.get('title', None)
    labels = kwargs.get('labels', ['Data', 'Cannon', 'BTSettl'])
    save   = kwargs.get('save', True)
    output = kwargs.get('out', 'Plot_Bands.pdf')

    if plot_band == 'OH':
        bands = [[15400,15450], [16350,16360], [16860,16890]]
    elif plot_band == 'Ca':
        bands = [[16131,16141], [16145,16155], [16152,16162]]
    elif plot_band == 'K':
        bands = [[15158,15168], [15163,15173]]
    elif plot_band == 'Mg':
        bands = [[15735,15745], [15743,15753], [15760, 15770]]
    elif plot_band == 'Al':
        bands = [[16713, 16723], [16745,16755], [16758,16768]]
    elif plot_band == 'Cannon':
        bands = [[15650,15780], [16150,16280]]
    elif plot_band == 'Full':
        bands = [[15200,15800],[15870,16420],[16490,16940]]
    #Regions from Rajpurohit paper:
    elif plot_band == 'R1':
        bands = [[15150,15450]]
    elif plot_band == 'R2':
        bands = [[15450,15800]]
    elif plot_band == 'R3':
        bands = [[15850,16420]]
    elif plot_band == 'R4':
        bands = [[16500,16910]]
    nbands = len(bands)

    fig, axs = plt.subplots(1, nbands, figsize=(16,4))

    nspecs = len(specs)
    colors = ['k', 'r', 'b', 'g']

    for i, ax in enumerate(fig.axes):
        for j in range(nspecs):
            ax.plot(specs[j].wave, specs[j].flux, color=colors[j], alpha=.8)
        if i==0:
            ax.set_ylabel(r'$F_{\lambda}$ [$erg/s \cdot cm^{2}$]', fontsize=20)
        ax.set_xlabel(r'$\lambda$ [$\mathring{A}$]', fontsize=20)
        ax.set_xlim(bands[i])
        ax.set_ylim([0.7, 1.15])
 
    if title != None:
        plt.suptitle(title, fontsize=25)

    if save == True:
        plt.savefig('plots/'+str(out))

    plt.show()
    plt.close()


def plotRajBands(**kwargs):

    #required
    specs = kwargs.get('specs')
    title = kwargs.get('title', 'Spectrum_Raj_Bands')
    
    #optional
    labels = kwargs.get('labels', ['Data', 'Cannon', 'BTSettl'])
    save   = kwargs.get('save', True)
    # output = kwargs.get('out', 'Plot_Bands.pdf')

    nspecs = len(specs)
    colors = ['k', 'r', 'b', 'g']
    bands  = [[15200,15450],[15450,15800],[15850,16420],[16500,16910]]

    fig, axs = plt.subplots(4, 1, figsize=(20,16), sharey=True) 

    for i, ax in enumerate(fig.axes):
        for j in range(nspecs):
            ax.plot(specs[j].wave, specs[j].flux, color=colors[j], alpha=.8, linewidth=1)
        ax.set_xlim(bands[i])
        if i == 0:
            ax.set_title(title, fontsize=20)

    if save == True:
        plt.savefig('plots/'+title+'.pdf')

    plt.show()
    plt.close()


# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# majorLocator = MultipleLocator(major)
#     majorFormatter = FormatStrFormatter('%d')
#     minorLocator = MultipleLocator(minor)

# ax.xaxis.set_major_locator(majorLocator)
#     ax.xaxis.set_major_formatter(majorFormatter)


