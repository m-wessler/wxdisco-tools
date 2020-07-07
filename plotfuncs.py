import gc
import numpy as np
import xarray as xr

from funcs import *
from config import *
from colortables import *

def mkdir_p(ipath):
    from os import makedirs, path
    import errno
    
    try:
        makedirs(ipath)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and path.isdir(ipath):
            pass
        else:
            raise
        
    return ipath

def xrsum(xarr):
    return xarr.sum(dim='time')

def load_mapdata(init_time, load_var, temp=True):
    from os import path, stat
    from multiprocessing import Pool, cpu_count

    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')
    
    dataset = []    
    for i, model in enumerate(models):
        for j in range(mcount): 
            # Add future support for other model downloads here
            if ensemble == 'sref':
                if j == 0:
                    member = model + '_ctl'
                elif j <= (mcount-1)/len(models):
                    member = model + '_n%i'%j
                elif j > (mcount-1)/len(models):
                    member = model + '_p%i'%(j-((mcount-1)/len(models)))

            elif ensemble == 'naefs':
                if j == 0:
                    member = model + '_c00'
                else:
                    member = model + '_p%02d'%j
            
            if temp == True:
                df = tmpdir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
                date, ensemble, dateH, dateH, member)
            else:
                df = datadir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
                    date, ensemble, dateH, dateH, member)
            
            dataset.append(
                xr.open_dataset(df, 
                decode_times=False,
                drop_variables=[v for v in output_vars 
                                if v not in load_var])[load_var])
        
    with Pool(cpu_count()-1) as p:
        acc_data = p.map(xrsum, dataset)
        p.close()
        p.join()

    [f.close() for f in dataset]
    [f.close() for f in acc_data]
    del dataset
    gc.collect()
    
    acc_data = xr.concat(acc_data, dim='member').load()    
    acc_data.rename('acc_' + load_var)
            
    return acc_data

def calc_stats(acc_data, init_time):
    dd = {}
    
    dd['init'] = init_time
    dd['lat'] = acc_data.lat
    dd['lon'] = acc_data.lon
    dd['mcount'] = acc_data.member.size

    dd['max'] = acc_data.max(dim='member')
    dd['min'] = acc_data.min(dim='member')
    dd['mean'] = acc_data.mean(dim='member')
    dd['med'] = acc_data.median(dim='member')
    dd['stdev'] = acc_data.std(dim='member')

    thresholds = snowthresh if 'snow' in acc_data.name else qpfthresh

    probs = [((xr.where(acc_data > thresh, 1, 0).sum(dim='member')/
          acc_data.member.size)*100) for thresh in thresholds]

    for k, v in zip(thresholds, probs):
        nk = 'prob_' + str(k)
        dd[nk] = v
        
    return dd

def load_plumedata(init_time, site, temp=False):
    from os import path, stat
    from multiprocessing import Pool, cpu_count

    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')

    print('Loading plume data for %s from file'%site)
    
    _dataset = []    
    for i, model in enumerate(models):
        for j in range(mcount): 
            # Add future support for other model downloads here
            if ensemble == 'sref':
                if j == 0:
                    member = model + '_ctl'
                elif j <= (mcount-1)/len(models):
                    member = model + '_n%i'%j
                elif j > (mcount-1)/len(models):
                    member = model + '_p%i'%(j-((mcount-1)/len(models)))

            elif ensemble == 'naefs':
                if j == 0:
                    member = model + '_c00'
                else:
                    member = model + '_p%02d'%j
            
            if temp == True:
                df = tmpdir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
                    date, ensemble, dateH, dateH, member)
            else:
                df = datadir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
                    date, ensemble, dateH, dateH, member)
                
            sitelat, sitelon = point_locs[site]

            _dataset.append(xr.open_dataset(df,
                decode_times=False))
            
    lons = _dataset[0].lon.values
    lats = _dataset[0].lat.values
    
    idx1d = (np.abs(lons-sitelon) + np.abs(lats-sitelat))
    idx = np.unravel_index(np.argmin(idx1d, axis=None), idx1d.shape)
    
    dataset = xr.concat(
        [mem.isel(x=idx[1], y=idx[0]).load() for mem in _dataset], 
        dim='member').compute()
    
    # Close the files immediately to free up the read for multiprocessing
    [f.close() for f in _dataset]

    dataset['acc_qpf'] = dataset.dqpf.cumsum(dim='time')
    dataset['acc_snow'] = dataset.snow.cumsum(dim='time')

    return dataset

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

def custom_violin(xr_data, axis, scale=1):
    vio_data = xr_data.values

    parts = axis.violinplot(
        vio_data, showmeans=False, showmedians=False,
        showextrema=False, widths=0.75)

    for ii, pc in enumerate(parts['bodies']):
        if ii == 0:
            iqrlab = 'IQR'
            medlab = 'Median'
        elif ii == 1:
            iqrlab, medlab = '', ''
        
        kdecolor = [0,0.4,0] if xr_data.name == 'dqpf' else [0,0,0.6]
        pc.set_facecolor(kdecolor)
        pc.set_edgecolor(kdecolor)
        pc.set_alpha(.5)

        quartile1, medians, quartile3 = np.percentile(
            vio_data, [25, 50, 75], axis=0)
        hiex, loex = np.percentile(vio_data, [5, 95], axis=0)
        whiskers = np.array([
            adjacent_values(sorted_array, q1, q3)
            for sorted_array, q1, q3 in zip(vio_data, quartile1, quartile3)])
        whiskersMin, whiskersMax = hiex, loex

        inds = np.arange(1, len(medians) + 1)
        axis.vlines(inds, quartile1, quartile3, color='k', linestyle='-', 
            lw=3*scale, label=iqrlab)
        axis.scatter(inds, medians, marker='_', color='r', s=25, zorder=3, 
            linewidth=2*scale, label=medlab)
        axis.vlines(inds, whiskersMin, whiskersMax, color='r', 
            linestyle='-', lw=1*scale, label='')
    
    axis.scatter(0,0, color='r', marker='|', label='5th/95th Percentile')
    axis.plot(0, 0, c=kdecolor, linewidth=5*scale, label='Density')

    return axis

def slr_spreadplot(data, axis, scale=1):
    slr = data.slr 

    slr_iqr = slr.quantile([0, .25, 0.50, .75, 1], 
        dim='member', interpolation='linear')

    axslr = axis.twinx()
    inds = range(1, slr.time.size+1)

    axslr.plot(inds, slr_iqr.sel(quantile=.50), 
                color='gray', alpha=0.75, zorder=-100, 
                linewidth=3.5*scale, label=('Median SLR'))

    # MEDIAN, IQR LINES
    # axslr.plot(inds, slr_iqr.sel(quantile=.75), 
    #             color='blue', alpha=0.75, zorder=-100, 
    #             linewidth=3.5*scale, label=('75th Percentile SLR'))

    # axslr.plot(inds, slr_iqr.sel(quantile=.25), 
    #             color='red', alpha=0.75, zorder=-100, 
    #             linewidth=3.5*scale, label=('25th Percentile SLR'))

    # MEDIAN, IQR SHADE
    axslr.fill_between(inds, 
        slr_iqr.sel(quantile=0.25), slr_iqr.sel(quantile=0.75),
        color='gray', alpha=0.4, zorder=-100, 
        linewidth=3*scale, label=('IQR'))

    # Median +/- STDEV
    # mslr = slr.median(dim='member')
    # slrstd = slr.std(dim='member')
    #     axslr.fill_between(inds, mslr - slrstd, mslr + slrstd,
    #         color='gray', alpha=0.4, zorder=-100, 
    #         linewidth=3*scale, label=(r'$\pm$ 1 SD'),)

    axslr.set_ylim(bottom=0,top=30)
    axslr.set_ylabel('Snow:Liqud Ratio')
    leg = axslr.legend(loc='upper right', fontsize='x-small')
   
    for line in leg.get_lines():
        line.set_linewidth(3.0)
   
    return None

def ensemble_plume(site, init_time, scale=0.75, temp=True):
    from os import remove as rm
    from os import system as sys

    import matplotlib
    matplotlib.use('agg')
    
    import matplotlib.pyplot as plt
    import matplotlib.dates as dates
    import matplotlib.image as mimage

    from pandas import to_datetime
    from subprocess import call as sys
    from datetime import datetime, timedelta

    data = load_plumedata(init_time, site, temp=temp)

    print('Plotting plume for %s'%site)
    
    sitelat, sitelon = point_locs[site]

    elev = xr.open_dataset(terrainfile)
    elats = elev.latitude.values
    elons = elev.longitude.values

    ilat = np.argmin(np.abs(elats - sitelat))
    ilon = np.argmin(np.abs(elons - sitelon))

    siteelev = elev.elevation.isel(latitude=ilat, longitude=ilon).values
    siteelev = int(round(siteelev * 3.28084, 0))
    elev.close()
    del elev

    plt.style.use('classic')
    matplotlib.rcParams.update({'font.size': 22*scale})

    fig, axs = plt.subplots(2, 2, figsize=(28*scale, 18*scale), 
        facecolor='w', edgecolor='k')

    fig.subplots_adjust(hspace=0.54, left=0.15, 
        bottom=0.13, right=0.96, top=0.88)

    axs = axs.flatten()

    if site in point_names.keys():
        plt.suptitle(('{} Downscaled Guidance at {:s} ({:s}: ' + 
                        '{:.2f} N {:.2f} W {:d} ft AMSL)\n' +
                        'Model Run: {} ({:d} Members)').format(
                        ensemble.upper(), point_names[site], 
                        site, sitelat, sitelon, siteelev,
                        init_time.strftime("%HZ %Y-%m-%d"), 
                        data.member.size, fontsize=14))

    else:
            plt.suptitle(('{} Downscaled Guidance at ' + 
                        '{:s} ({:.2f} N {:.2f} W {:d} ft AMSL)\n' + 
                        'Model Run: {} ({:d} Members)').format(
                        ensemble.upper(), site, sitelat, sitelon, siteelev,
                        init_time.strftime("%HZ %Y-%m-%d"), 
                        data.member.size, fontsize=14))

    half = int(data.member.size / 2)

    dtarr = to_datetime(data.time.values)

    for i, member in data.groupby('member'):
        lw, mlw = 1.5*scale, 5*scale

        if member.member == 1:
            label = (str(member.member_id.values).split('_')[0].upper() + 
                ' Members')
            qcol = [0,0.9,0]
            scol = [0,0,1.0]

        elif member.member == half+1:
            label = (str(member.member_id.values).split('_')[0].upper() + 
                ' Members')
            qcol = [0.7,0.9,0]
            scol = [0,0.6,0.9]

        else:
            label = ''

        axs[0].plot(dtarr, member.acc_qpf,
                    color=qcol, linewidth=lw, label=label)
        axs[2].plot(dtarr, member.acc_snow,
                    color=scol, linewidth=lw, label=label)

    halflab0 = str(data.member_id[:half][0].values).split('_')[0].upper()
    halflab1 = str(data.member_id[half:][0].values).split('_')[0].upper()

    axs[0].plot(dtarr, data.acc_qpf[:half].mean(dim='member'),
                '-', color='r', linewidth=mlw, label='%s Mean'%halflab0)
    axs[0].plot(dtarr, data.acc_qpf[half:].mean(dim='member'),
                '--', color='r', linewidth=mlw, label='%s Mean'%halflab1)
    axs[0].plot(dtarr, data.acc_qpf[:].mean(dim='member'),
                '-', color=[0,0.4,0], linewidth=mlw, 
                label='%s Mean'%ensemble.upper())

    axs[2].plot(dtarr, data.acc_snow[:half].mean(dim='member'),
                '-', color='r', linewidth=mlw, label='%s Mean'%halflab0)
    axs[2].plot(dtarr, data.acc_snow[half:].mean(dim='member'),
                '--', color='r', linewidth=mlw, label='%s Mean'%halflab1)
    axs[2].plot(dtarr, data.acc_snow[:].mean(dim='member'),
                '-', color=[0,0,0.6], linewidth=mlw,
                label='%s Mean'%ensemble.upper())

    axs[1] = custom_violin(data.dqpf, axs[1], scale=scale)
    axs[3] = custom_violin(data.snow, axs[3], scale=scale)

    slr_spreadplot(data, axs[3], scale=scale)

    subtitles = ['Accumulated Precipitation',
                '%d-Hourly Precipitation'%fhrstep,
                'Accumulated Snow',
                '%d-Hourly Snow'%fhrstep]

    ylabels = ['Precipitation (Inches)', 'Precipitation (Inches)', 
                'Snow (Inches)', 'Snow (Inches)']

    logo = mimage.imread(scriptdir + '../Ulogo_400p.png')
    axs_water = [[0.16, 0.73, .065, .11],
                 [0.89, 0.73, .065, .11],
                 [0.16, 0.27, .065, .11],
                 [0.89, 0.27, .065, .11]]
    
    for i, ax in enumerate(axs):
        ax.set_title(subtitles[i])
        ax.set_xlabel('Day/Hour')
        ax.set_ylabel(ylabels[i])

        if i%2 == 0:
            _ytop = 0.5 if i == 0 else 1.0
            ytop = (ax.get_ylim()[1] * 1.2 
                if ax.get_ylim()[1] >= _ytop else _ytop)

            ax.set_ylim([0, ytop])
            ax.set_xlim([dtarr.min(), dtarr.max()])
            
            if ensemble == 'sref':
                ax.xaxis.set_major_locator(dates.HourLocator(
                    byhour=range(0, 24, 6)))
                ax.xaxis.set_major_formatter(dates.DateFormatter('%d/%HZ'))
                plt.setp(ax.xaxis.get_majorticklabels(), rotation=80)

            leg = ax.legend(loc=2, handlelength=5*scale, 
                fontsize='x-small', ncol=3)
            ax.grid(b=True, linestyle='dotted')

        else:
            ytop = ax.get_ylim()[1] * 1.05 if ax.get_ylim()[1] >= 0.5 else 0.5
            ax.set_ylim([0, ytop])

            tick_hack = [(dti.strftime('%d/%HZ') if dti.hour % 6 == 0
                        else '') for dti in dtarr]
            
            if ensemble == 'sref':
                ax.set_xticks(np.arange(1, len(tick_hack)+1))
                ax.set_xticklabels(tick_hack, rotation=80)
                
            ax.set_xlim([1.5, len(tick_hack) + 0.5])

            leg = ax.legend(loc=2, handlelength=5*scale, markerscale=5*scale, 
                scatterpoints=1, fontsize='x-small', ncol=2)
            ax.grid(b=True, axis='y', linestyle='dotted')
        
        if ensemble == 'naefs':
            disptype = 'date'
            monthlist = ['Jan','Feb','Mar','Apr','May','Jun',
                            'Jul','Aug','Sep','Oct','Nov','Dec']

            ticklabs = []
            ticksites = np.arange(fhrstart, fhrend + fhrstep, fhrstep)

            for t in ticksites:
                d = init_time + timedelta(float(t) / 24.0)
                ticklab = ''
                
                if np.floor(float(t)/12.0) == float(t) / 12.0:
                    ticklab = '%02d'%d.hour + 'z'
                    
                if np.floor(float(t)/24.0) == float(t) / 24.0:
                    ticklab = (ticklab + '\n%02d'%d.day + '-' + 
                        monthlist[d.month - 1])
                        
                ticklabs.append(ticklab)

            ticklocs = dtarr if i%2 == 0 else np.arange(1, len(tick_hack) + 1)
            ax.set_xticks(ticklocs)
            ax.set_xticklabels(ticklabs)

        for line in leg.get_lines():
            line.set_linewidth(8*scale)

        ax_water = fig.add_axes(axs_water[i])
        ax_water.axis('off')
        ax_water.imshow(logo, aspect='auto', zorder=10, alpha=0.35)

    date = init_time.strftime('%Y%m%d')
    dateH = init_time.strftime('%Y%m%d%H')
    figdir = mkdir_p(imgdir + '%s/images/models/%s/'%(date, ensemble))

    figpath = figdir + '{}PL_{}{}F{:03d}.png'.format(
        ensemble.upper(), site, dateH, int(fhrend-fhrstart))
    figpath_gif = figpath[:-4] + '.gif'

    plt.savefig(figpath, bbox_inches='tight')
    plt.close()

    sys('convert ' + figpath + ' ' + figpath_gif, shell=True)
    rm(figpath)

    return None

def ensemble_qpf6panel(region, dd):
    import os
    from scipy import ndimage
    from os import remove as rm
    from subprocess import call as sys
    from datetime import datetime, timedelta
    
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mimage
    from mpl_toolkits.basemap import Basemap, maskoceans
    
    # Other very nice styles available. See the docs
    plt.style.use('classic')

    print('Plotting QPF 6-panel for %s'%region)

    # Gather the metadata needed
    init_time = dd['init']
    fx_len = int(fhrend)-int(fhrstart)
    
    # Crop the desired region
    lons, lats = dd['lon'].values, dd['lat'].values
    minlon, maxlon, minlat, maxlat = map_regions[region]

    # Initialize the figure frame
    fig = plt.figure(num=None, figsize=(11.0*scale,8.0*scale), dpi=300*scale, 
                     facecolor='w', edgecolor='k')

    # Customize the subplots
    fig.subplots_adjust(top=0.89, bottom=0.02, left=0.055, 
                        right=0.99, hspace=0.1*scale, wspace=0.22)

    plt.suptitle(('{} Downscaled QPF ({:d} members) Initialized {}\n' + 
        '{:d}-{:d}-h Forecast Valid {} to {}').format(
            ensemble.upper(),
            dd['mcount'], 
            dd['init'].strftime('%HZ %Y-%m-%d'), 
            fhrstart, fhrend, 
            dd['init'].strftime('%HZ %Y-%m-%d'), 
            (dd['init'] + timedelta(
                hours=fx_len)).strftime(
                '%HZ %Y-%m-%d')), fontsize=14)
    
    # Initialize the basemap object
    bmap = Basemap(
        llcrnrlon=minlon-0.01,
        llcrnrlat=minlat-0.01,
        urcrnrlon=maxlon+0.01,
        urcrnrlat=maxlat+0.01,
        resolution='i',
        # Project as 'aea', 'lcc', or 'mill'
        projection='mill')

    # Project the lat/lon grid
    x, y = bmap(lons, lats)

    # Types of plots to draw
    plottype = (['mean', 'max', 'min'] + 
                ['prob_' + str(threshold) 
                 for threshold in qpfthresh])

    # Titles for the plots (indexes must align with above)
    titles = ['Ensemble Mean %d-hr Accum Precip'%fx_len,
                'Ensemble Max %d-hr Accum Precip'%fx_len,
                'Ensemble Min %d-hr Accum Precip'%fx_len,
                'Prob of %d-hr Precip > %s" (%%)'%(fx_len, qpfthresh[0]),
                'Prob of %d-hr Precip > %s" (%%)'%(fx_len, qpfthresh[1]),
                'Prob of %d-hr Precip > %s" (%%)'%(fx_len, qpfthresh[2])]

    # Painstakingly placed colorbars...
    cbar_axes = [[0.01,0.52,0.012,0.32],
                 [0.34,0.52,0.012,0.32],
                 [0.67,0.52,0.012,0.32],
                 [0.01,0.08,0.012,0.32],
                 [0.34,0.08,0.012,0.32],
                 [0.67,0.08,0.012,0.32]]

    # Painstakingly placed logo...
    logo = mimage.imread(scriptdir + '../Ulogo_400p.png')
    bot, top, left = .05, 0.505, 0.245
    axs_water = [[left, top, .065, .11],
                 [left+0.332, top, .065, .11], 
                 [left+0.664, top, .065, .11], 
                 
                 [left, bot, .065, .11], 
                 [left+0.332, bot, .065, .11],
                 [left+0.664, bot, .065, .11]]
    
    # If choosing to plot elevation contours
    if plot_elevation == True:
        zsmooth, zint = elev_smoothing[region]
        
        # Contour interval
        zlevs = np.arange(0, 4500+.01, zint)
        
        # Smoothing parameters
        efold = (res_ensemble/zsmooth) / res_prism + 1
        sigma = efold / (np.pi*np.sqrt(2))
        
        # Smooth the terrain DEM as desired
        z = ndimage.filters.gaussian_filter(
            get_elev(dd[plottype[0]]), 
            sigma, mode='nearest')
    
    # List to store the colorbar data
    cbd = []
    
    for i in range(6):
        # Basemap likes explicit subplots
        ax = plt.subplot(2,3,i+1)
        ax.set_title(titles[i], fontsize=11)

        bmap.drawcoastlines(linewidth=1.0, color='black')
        bmap.drawcountries(linewidth=0.85, color='black')
        bmap.drawstates(linewidth=1.25, color='black')
        
        if plot_elevation == True:
            bmap.contour(x, y, z, levels=zlevs, colors='black', alpha=0.75, 
            linewidths=1)

        # Choose the appropriate colormaps and levels for the plot
        # Edit these in colortables.py
        if 'prob' in plottype[i]:
            levs, norml, cmap, tickloc, ticks = (problevs, probnorml, 
                                                 probcmap, problevloc, 
                                                 probticks)
        else:
            levs, norml, cmap, tickloc, ticks = (qpflevs, qpfnorml, 
                                                 qpfcmap, qpflevloc, qpfticks)

        # Mask the oceans, inlands=False to prevent masking lakes (e.g. GSL)
        data_lsmask = maskoceans(lons, lats, dd[plottype[i]], inlands=False)

        # Allow for precip beyond the maximum contourf
        extend = 'max' if 'prob' not in plottype[i] else 'neither'
        cbd.append(bmap.contourf(x, y, data_lsmask, extend=extend,
                                 levels=levs, norm=norml, cmap=cmap, 
                                 alpha=0.80))

        # Manually establish the colorbar
        cax = fig.add_axes(cbar_axes[i])
        cbar = plt.colorbar(cbd[i], cax=cax, orientation='vertical') 
        cbar.set_ticks(tickloc)
        cbar.set_ticklabels(ticks)
        cbar.ax.tick_params(size=0)
        cbytick_obj = plt.getp(cbar.ax.axes, 'yticklabels')
        plt.setp(cbytick_obj, fontsize=8*scale, rotation=0)

        # Add the U watermark
        ax_water = fig.add_axes(axs_water[i])
        ax_water.axis('off')
        ax_water.imshow(logo, aspect='auto', zorder=10, alpha=0.40)
  
        # Determine the max precip in the visible frame
        ld = dd[plottype[i]]
        if 'prob' not in plottype[i]:
            ld_mask = xr.where(
                ((minlat <= ld.lat) & (ld.lat <= maxlat)) & 
                ((minlon <= ld.lon) & (ld.lon <= maxlon)), ld, np.nan)
            localmax = ld_mask.max()
            xlon, xlat = np.where(ld_mask == localmax)
              
            # Mark the max precip in the frame with an 'x' (Off for now)
            # Basemap is improperly placing the lat lon point... why?
            # xx, yy = x[xlon, xlat], y[xlon, xlat]
            # bmap.scatter(xx, yy, marker='x', c='k', s=100)
            
            ax.set_xlabel('Max: {:.02f} in ({:.02f}, {:.02f})'.format(
                localmax.values, 
                lats[xlon, xlat][0], lons[xlon, xlat][0]), fontsize=8*scale)
            
        elif ((i == 4) & (plot_elevation == True)):
            # Annotate the bottom of the plot to explain the elevation contour
            ax.set_xlabel('Elevation (Gray Contours) Every 1000 m', 
            fontsize=8*scale)
    
    # Save the file according to previous convention
    date = init_time.strftime('%Y%m%d')
    dateH = init_time.strftime('%Y%m%d%H')
    figdir = mkdir_p(imgdir + '%s/images/models/%s/'%(date, ensemble))
        
    figpath = figdir + '{}PP_{}{}F{:03d}.png'.format(
        ensemble.upper(), region, dateH, int(fhrend-fhrstart))

    figpath_gif = figpath[:-4] + '.gif'

    plt.savefig(figpath, bbox_inches='tight')
    plt.close()
    
    # Convert to gif using imagemagick
    sys('convert ' + figpath + ' ' + figpath_gif, shell=True)
    rm(figpath)

    return None

def ensemble_snow6panel(region, dd):
    import os
    from scipy import ndimage
    from os import remove as rm
    from subprocess import call as sys
    from datetime import datetime, timedelta
    
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mimage
    from mpl_toolkits.basemap import Basemap, maskoceans
    
    # Other very nice styles available. See the docs
    plt.style.use('classic')

    print('Plotting Snow 6-panel for %s'%region)

    # Gather the metadata needed
    init_time = dd['init']
    fx_len = int(fhrend)-int(fhrstart)
    
    # Crop the desired region
    lons, lats = dd['lon'].values, dd['lat'].values
    minlon, maxlon, minlat, maxlat = map_regions[region]

    # Initialize the figure frame
    fig = plt.figure(num=None, figsize=(11.0*scale,8.0*scale), dpi=300*scale, 
                     facecolor='w', edgecolor='k')

    # Customize the subplots
    fig.subplots_adjust(top=0.89, bottom=0.02, left=0.055, 
                        right=0.99, hspace=0.1*scale, wspace=0.22)
    
    plt.suptitle(('{} Downscaled QPF/Snow ({:d} members) Initialized {}\n' + 
        '{:d}-{:d}-h Forecast Valid {} to {}').format(
            ensemble.upper(),
            dd['mcount'], 
            dd['init'].strftime('%HZ %Y-%m-%d'), 
            fhrstart, fhrend, 
            dd['init'].strftime('%HZ %Y-%m-%d'), 
            (dd['init'] + timedelta(
                hours=fx_len)).strftime(
                '%HZ %Y-%m-%d')), fontsize=14)

    # Initialize the basemap object
    bmap = Basemap(
        llcrnrlon=minlon-0.01,
        llcrnrlat=minlat-0.01,
        urcrnrlon=maxlon+0.01,
        urcrnrlat=maxlat+0.01,
        resolution='i',
        # Project as 'aea', 'lcc', or 'mill'
        projection='mill')

    # Project the lat/lon grid
    x, y = bmap(lons, lats)

    # Types of plots to draw
    plottype = (['mqpf', 'mean'] + 
                ['prob_' + str(threshold) 
                 for threshold in snowthresh])
    
    # Titles for the plots (indexes must align with above)
    titles = ['Ensemble Mean %d-hr Accum Precip'%fx_len,
                'Ensemble Mean %d-hr Accum Snow'%fx_len,
                'Prob of %d-hr Snow > %s" (%%)'%(fx_len, snowthresh[0]),
                'Prob of %d-hr Snow > %s" (%%)'%(fx_len, snowthresh[1]),
                'Prob of %d-hr Snow > %s" (%%)'%(fx_len, snowthresh[2]),
                'Prob of %d-hr Snow > %s" (%%)'%(fx_len, snowthresh[3])]

    # Painstakingly placed colorbars...
    cbar_axes = [[0.01,0.52,0.012,0.32],
                 [0.34,0.52,0.012,0.32],
                 [0.67,0.52,0.012,0.32],
                 [0.01,0.08,0.012,0.32],
                 [0.34,0.08,0.012,0.32],
                 [0.67,0.08,0.012,0.32]]

    # Painstakingly placed logo...
    logo = mimage.imread(scriptdir + '../Ulogo_400p.png')
    bot, top, left = .05, 0.505, 0.245
    axs_water = [[left, top, .065, .11],
                 [left+0.332, top, .065, .11], 
                 [left+0.664, top, .065, .11], 
                 
                 [left, bot, .065, .11], 
                 [left+0.332, bot, .065, .11],
                 [left+0.664, bot, .065, .11]]
    
    # If choosing to plot elevation contours
    if plot_elevation == True:
        zsmooth, zint = elev_smoothing[region]
        
        # Contour interval
        zlevs = np.arange(0, 4500+.01, zint)
        
        # Smoothing parameters
        efold = (res_ensemble/zsmooth) / res_prism + 1
        sigma = efold / (np.pi*np.sqrt(2))
        
        # Smooth the terrain DEM as desired
        z = ndimage.filters.gaussian_filter(
            get_elev(dd[plottype[0]]), 
            sigma, mode='nearest')
    
    # List to store the colorbar data
    cbd = []
    
    for i in range(6):
        # Basemap likes explicit subplots
        ax = plt.subplot(2,3,i+1)
        ax.set_title(titles[i], fontsize=11)

        bmap.drawcoastlines(linewidth=1.0, color='black')
        bmap.drawcountries(linewidth=0.85, color='black')
        bmap.drawstates(linewidth=1.25, color='black')
        
        # Choose the appropriate colormaps and levels for the plot
        # Edit these in colortables.py
        if plot_elevation == True:
            bmap.contour(x, y, z, 
            levels=zlevs, colors='black', alpha=0.75, linewidths=1)

        if 'prob' in plottype[i]:
            levs, norml, cmap, tickloc, ticks = (problevs, probnorml, 
                                                 probcmap, problevloc, 
                                                 probticks)
        else:
            if 'qpf' in plottype[i]:
                levs, norml, cmap, tickloc, ticks = (qpflevs, qpfnorml, 
                                                 qpfcmap, qpflevloc, qpfticks)
            else:
                levs, norml, cmap, tickloc, ticks = (snowlevs, snownorml, 
                                                     snowcmap, snowlevloc, 
                                                     snowticks)

        # Mask the oceans, inlands=False to prevent masking lakes (e.g. GSL)
        data_lsmask = maskoceans(lons, lats, dd[plottype[i]], inlands=False)

        # Allow for precip beyond the maximum contourf
        extend = 'max' if 'prob' not in plottype[i] else 'neither'
        cbd.append(bmap.contourf(x, y, data_lsmask, extend=extend,
                                 levels=levs, norm=norml, 
                                 cmap=cmap, alpha=0.90))

        # Manually establish the colorbar
        cax = fig.add_axes(cbar_axes[i])
        cbar = plt.colorbar(cbd[i], cax=cax, orientation='vertical') 
        cbar.set_ticks(tickloc)
        cbar.set_ticklabels(ticks)
        cbar.ax.tick_params(size=0)
        cbytick_obj = plt.getp(cbar.ax.axes, 'yticklabels')
        plt.setp(cbytick_obj, fontsize=8*scale, rotation=0)

        # Add the U watermark
        ax_water = fig.add_axes(axs_water[i])
        ax_water.axis('off')
        ax_water.imshow(logo, aspect='auto', zorder=10, alpha=0.40)
        
        # Determine the max precip in the visible frame
        ld = dd[plottype[i]]
        if 'prob' not in plottype[i]:
            ld_mask = xr.where(
                ((minlat <= ld.lat) & (ld.lat <= maxlat)) & 
                ((minlon <= ld.lon) & (ld.lon <= maxlon)), ld, np.nan)
            localmax = ld_mask.max()
            xlon, xlat = np.where(ld_mask == localmax)
              
            # Mark the max precip in the frame with an 'x' (Off for now)
            # Basemap is improperly placing the lat lon point... why?
            # xx, yy = x[xlon, xlat], y[xlon, xlat]
            # bmap.scatter(xx, yy, marker='x', c='k', s=100)
            
            ax.set_xlabel('Max: {:.02f} in ({:.02f}, {:.02f})'.format(
                localmax.values, 
                lats[xlon, xlat][0], lons[xlon, xlat][0]), fontsize=8*scale)
            
        elif ((i == 4) & (plot_elevation == True)):
            # Annotate the bottom of the plot to explain the elevation contour
            ax.set_xlabel('Elevation (Gray Contours) Every 1000 m', 
            fontsize=8*scale)
    
    # Save the file according to previous convention
    date = init_time.strftime('%Y%m%d')
    dateH = init_time.strftime('%Y%m%d%H')
    figdir = mkdir_p(imgdir + '%s/images/models/%s/'%(date, ensemble))

    figpath = figdir + '{}PS_{}{}F{:03d}.png'.format(
        ensemble.upper(), region, dateH, int(fhrend-fhrstart))
        
    figpath_gif = figpath[:-4] + '.gif'

    plt.savefig(figpath, bbox_inches='tight')
    plt.close()
    
    # Convert to gif using imagemagick
    sys('convert ' + figpath + ' ' + figpath_gif, shell=True)
    rm(figpath)

    return None