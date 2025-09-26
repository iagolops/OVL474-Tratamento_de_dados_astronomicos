import matplotlib.pyplot as plt
import numpy as np
#plt.style.use('/pscratch/sd/i/iago/iago.mplstyle')
colors = ['#5495CFFF', '#9E8ABCFF', '#5FB785FF', '#7C873EFF', '#F5AF4DFF', '#DB4743FF']

def plot_mag_hist_comparison(cat1, cat2,
                             title='Magnitude histograms (S/N > 10)',
                             bins=100, xlim=[18, 30],
                             sn_cut=10, mag_i_limit=24.5,
                             colors=('C0','C1'),  # e.g. ('C5','C4')
                             save=False, path=''):
    """
    Compare magnitude histograms (S/N > sn_cut) for two catalogs across LSST bands.

    Parameters
    ----------
    cat1, cat2 : DataFrame
        Must contain 'mag_{band}_lsst' and 'mag_err_{band}_lsst' for some bands.
    title : str
        Figure title.
    bins : int
        Number of bins.
    xlim : list[float, float]
        Magnitude axis limits.
    sn_cut : float
        Minimum S/N (1.086 / mag_err) to include an object.
    mag_i_limit : float
        Maximum i-band magnitude for preselection.
    colors : tuple
        Two colors for cat1 and cat2 histograms.
    save : bool
        Whether to save the figure.
    path : str
        File path to save.
    """
    bands = ['u','g','r','i','z','y']
    fig, axes = plt.subplots(3, 2, figsize=(12, 16))
    fig.suptitle(title, fontsize=16)

    cat1 = cat1.copy()
    cat2 = cat2.copy()
    
    sn_i_1 = 1.086 / cat1['mag_err_i_lsst']
    mask_1 = (sn_i_1 >= sn_cut) & (cat1['mag_i_lsst'] <= mag_i_limit) & (cat1['mag_i_lsst'] > 17.5)
    cat1 = cat1.loc[mask_1]

    sn_i_2 = 1.086 / cat2['mag_err_i_lsst']
    mask_2 = (sn_i_2 >= sn_cut) & (cat2['mag_i_lsst'] <= mag_i_limit) & (cat2['mag_i_lsst'] > 17.5)
    cat2 = cat2.loc[mask_2]
    
    for ax, band in zip(axes.flatten(), bands):
        mag_col = f'mag_{band}_lsst'
        err_col = f'mag_err_{band}_lsst'

        def get_mag(cat):
            if mag_col not in cat.columns or err_col not in cat.columns:
                return None
            sn = 1.086 / cat[err_col]
            mask = (sn > sn_cut) & (cat[mag_col] > xlim[0]) & (cat[mag_col] < xlim[1])
            return cat.loc[mask, mag_col].values

        mag1 = get_mag(cat1)
        mag2 = get_mag(cat2)

        if mag1 is not None and mag1.size > 0:
            ax.hist(mag1, bins=bins, range=xlim,
                    histtype='stepfilled', alpha=0.3,
                    color=colors[0], density=True, label='Roman Rubin (HSC error)')
            ax.hist(mag1, bins=bins, range=xlim,
                    histtype='step', color=colors[0],
                    linewidth=1.5, density=True)
        if mag2 is not None and mag2.size > 0:
            ax.hist(mag2, bins=bins, range=xlim,
                    histtype='stepfilled', alpha=0.3,
                    color=colors[1], density=True, label='HSC Y3 catalog')
            ax.hist(mag2, bins=bins, range=xlim,
                    histtype='step', color=colors[1],
                    linewidth=1.5, density=True)

        ax.set_xlabel(f"mag {band}", fontsize=14)
        ax.set_ylabel("Normalized counts", fontsize=14)
        ax.set_xlim(xlim)
        ax.set_ylim(1e-3, 1e0)
        ax.set_yscale('log')
        ax.grid(True, linestyle='--', alpha=0.6)
        if band == 'g':
            ax.legend(fontsize=12, loc='center right')

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save and path:
        plt.savefig(path, dpi=150)
    plt.show()


def plot_errors(catalog, title='Errors', gridsize=[400, 200], bins='log', cmap='inferno',
                xlim=[20, 30], ylim=[0, 100], save=False, path=''):
    """
    Plots the signal-to-noise ratio (S/N) as a function of magnitude for HSC like bands (g,r,i,z,y).

    Parameters:
    - catalog: dictionary or DataFrame with 'mag_{band}_lsst' and 'mag_err_{band}_lsst' columns
    - title: main plot title
    - gridsize: resolution of the hexbin plot
    - bins: scale for the hexbin color mapping ('log', 'linear', etc.)
    - cmap: colormap used in hexbin
    - xlim, ylim: axis limits
    - save: whether to save the figure
    - path: path to save the figure if save is True
    """

    bands = ['g', 'r', 'i', 'z', 'y'] 
    fig, axes = plt.subplots(3, 2, figsize=[12, 16])
    fig.suptitle(title, fontsize=16)

    catalog = catalog[catalog['mag_i_lsst']<24.5]

    i = 0

    for ax, band in zip(axes.flatten(), bands):
        mag = np.array(catalog[f'mag_{band}_lsst'])
        err = np.array(catalog[f'mag_err_{band}_lsst'])

        sn = 1.086/err
        
        mask = (sn < 100) & (mag < 30) & (mag > 20)
        catalog = catalog[mask]

        # Recompute magnitude and S/N after applying the mask
        mag = np.array(catalog[f'mag_{band}_lsst'])
        err = np.array(catalog[f'mag_err_{band}_lsst'])
        sn = 1.086/err

        ax.hexbin(mag, sn, gridsize=gridsize, cmap=cmap, bins=bins, mincnt=1)
        ax.set_ylabel("S/N", fontsize=18)
        ax.set_xlabel(f"mag {band}", fontsize=18)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        
        ax.axhline(10, color='red', label='10σ', linestyle='--')

        ax.grid(True, linestyle='--', alpha=0.6)
        ax.legend(fontsize=18, loc='upper right')
        ax.tick_params(axis='both', which='major', labelsize=12)

    fig.tight_layout(rect=[0, 0, 1, 0.97])

    if save and path:
        plt.savefig(path, dpi=150)
    plt.show()


def plot_diff_errors(cat1, cat2, gridsize=[100,100], 
                     xlim=[20, 28], ylim=[0, 100], years=[], save=False, path=''):
    bands = ['u', 'g', 'r', 'i', 'z', 'y']

    fig, axes = plt.subplots(2, 3, figsize=[18, 10], sharey=True)
    im = None
    for i, band in enumerate(bands):
        ax = axes.flatten()[i]

        def compute_grid(catalog):
            if f'mag_{band}_lsst' not in catalog or f'mag_err_{band}_lsst' not in catalog:
                return None

            mask = (catalog['mag_i_lsst'] < 24.5) & (1.086 / catalog['mag_err_i_lsst']>10)
            catalog = catalog[mask]

            mag = np.array(catalog[f'mag_{band}_lsst'])
            err = np.array(catalog[f'mag_err_{band}_lsst'])
            sn = 1.086 / err

            mask = (sn < ylim[1]) & (sn > ylim[0]) & (mag < xlim[1]) & (mag > xlim[0])
            mag = mag[mask]
            sn = sn[mask]

            return np.histogram2d(mag, sn, bins=gridsize, range=[xlim, ylim], density=True)[0]

        h1 = compute_grid(cat1)
        h2 = compute_grid(cat2)
        if h1 is None:
            ax.axis('off')
            continue
        if h2 is None:
            h2 = np.zeros_like(h1)

        diff = (h1 - h2)
        extent = [xlim[0], xlim[1], ylim[0], ylim[1]]
        im = ax.imshow(diff.T, extent=extent, origin='lower', aspect='auto',
                       cmap='coolwarm', vmin=-0.05, vmax=0.05)

        ax.set_xlabel(f'{band}-band', fontsize=17)
        if band=='u' or band=='i':
            ax.set_ylabel('S/N', fontsize=18)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xticks([21, 23, 25, 27])      
        ax.set_xticklabels(['21', '23', '25', '27'])
        ax.text(23.5,80,f'nYrObs: {years[i]}', fontsize=18)
        ax.tick_params(top=False, right=False)
        
        col = i % 3
        if col > 0:
            ax.set_yticklabels([])

    fig.subplots_adjust(wspace=0, hspace=0.3)

    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.97)
    cbar.set_label(r'$N_{\rm rubin} - N_{\rm HSC}\,/\,N_{\rm tot}$', fontsize=20)


    if save and path:
        plt.savefig(path, dpi=500)
    plt.show()