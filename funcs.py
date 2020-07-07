import gc
import numpy as np
import xarray as xr

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

def bytes2human(n):
    ''' http://code.activestate.com/recipes/578019 '''

    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n

def get_init(req_time=None):
    from sys import exit
    from datetime import datetime, timedelta
    
    if req_time is not None:
        try:
            mostrecent = datetime.strptime(req_time, '%Y%m%d%H')
        except:
            print('Invalid time requested, please enter as YYYYMMDDHH')
            exit()

    else:
        qinit = datetime.utcnow() - timedelta(hours=delay)
        qinitH = qinit.hour

        inits = naefs_inits if ensemble == 'naefs' else sref_inits
        
        mostrecent = None
        for i in inits:
            if qinitH >= i:
                mostrecent = qinit.strftime('%Y%m%d') + str(i)
                mostrecent = datetime.strptime(mostrecent, '%Y%m%d%H')
                break

        if mostrecent is None:
            # When hour outside of the above cases, must catch exception
            mostrecent = ((qinit - timedelta(days=1)).strftime('%Y%m%d') 
                + str(inits[0]))
            mostrecent = datetime.strptime(mostrecent, '%Y%m%d%H')
    
    return mostrecent

def get_grids(init_time):    
    from os import path, stat, remove
    from multiprocessing import Pool, cpu_count, get_context
    from datetime import datetime

    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')
    
    gribdir = mkdir_p(tmpdir + '%s/models/%s/'%(date, ensemble))
    
    dl_list = []
    dl_start = datetime.utcnow()
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
            
            for fhr in range(fhrstart, fhrend+1, fhrstep):
                uuname = '%s%sF%02i.grib2'%(dateH, member, fhr)
                checkfile = gribdir + uuname
                
                # Check if file exists
                if (path.isfile(checkfile)):
                    pass

                    # Check for valid filesize (junk data)
                    if (stat(checkfile).st_size < minsize):
                        # If not, download grids to disk
                        dl_list.append([init_time, fhr, member, checkfile])
                        remove(checkfile)
                    else:
                        pass
                
                else:
                    # If not, download grids to disk
                    dl_list.append(
                        [init_time, fhr, member, checkfile, dl_start])

    # Since most of the download verification is now happening within
    # the worker pool, this isn't actually being utilized as a while loop
    # should right now. Can either finish setting this up to only delete a 
    # sucessful return from the list or remove the loop entirely...
    # though leaving as-is is not a problem either.
    while len(dl_list) > 0:
        cores = len(dl_list) if len(dl_list) < cpu_count()-1 else cpu_count()-1
        
        if mpi_limit is not None:
            cores = cores if cores <= mpi_limit else mpi_limit
        
        print('Downloading %i files on %i cores'%(len(dl_list), cores))
        
        with get_context(spawntype).Pool(cores) as p:
            post_dl_list = p.map(download_grib, dl_list)
            p.close()
            p.join()
        
        [dl_list.remove(dl) for dl in post_dl_list]
        del post_dl_list

    print('Found all files for %s %s'%(ensemble, init_time))
                
    return None

def download_grib(params):
    # Download the desired model run and verify
    import urllib.request
    from sys import exit
    from time import sleep
    from datetime import datetime
    from os import remove, stat, path
    from subprocess import call

    date = params[0].strftime('%Y%m%d')
    ihr = params[0].hour
    fhr = params[1]
    family = params[2].split('_')[0]
    member = params[2].split('_')[1]
    fpath = params[3]
    start_time = params[4]
    
    # Add future support for other model downloads here
    # Use NCDC NOMADS filterscripts
    # http://nomads.ncep.noaa.gov/
    if ensemble == 'sref':
        base = ('https://nomads.ncep.noaa.gov/cgi-bin/' +
            'filter_sref_132.pl?file=sref')
        mid = '_{}.t{:02d}z.pgrb132.{}.f{:02d}.grib2'.format(
            family, ihr, member, fhr)
        webdir = '&dir=%2Fsref.{}%2F{:02d}%2Fpgrb'.format(date, ihr)
                
    elif ensemble == 'naefs':
        if family == 'gefs':

            base = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gens_0p50.pl'
            mid = '?file=ge{}.t{:02d}z.pgrb2a.0p50.f{:03d}'.format(
                member, ihr, fhr)
            webdir = '&dir=%2Fgefs.{}%2F{:02d}%2Fpgrb2ap5'.format(date, ihr)

        elif family == 'cmce':
            base = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_cmcens.pl'
            mid = '?file=cmc_ge{}.t{:02d}z.pgrb2a.0p50.f{:03d}'.format(
                member, ihr, fhr)
            webdir = '&dir=%2Fcmce.{}%2F{:02d}%2Fpgrb2ap5'.format(date, ihr)
    
    mvars = '&var_APCP=on&var_HGT=on&var_TMP=on&var_RH=on'

    mlevs = ('&lev_500_mb=on&lev_700_mb=on&lev_850_mb=on' + 
    '&lev_925_mb=on&lev_1000_mb=on&lev_surface=on')

    subset = ('&subregion=&leftlon={}&rightlon={}'.format(
        minlon, maxlon) + '&toplat={}&bottomlat={}'.format(
        maxlat, minlat))

    url = base + mid + mvars + mlevs + subset + webdir

    # Download the grib to disk
    while not path.isfile(fpath):
        try:
            urllib.request.urlretrieve(url, fpath)

        except OSError:
            # Sometimes urllib struggles. Before totally giving up, try this
            # the old fashioned way first...
            curlcommand = 'curl -s -m {} -o {} {}'.format(timeout, fpath, url)
            call(curlcommand, shell=True)

        try:
            fsize = stat(fpath).st_size
        except:
            print('FILE NOT FOUND Data not yet available. Waiting', 
                wait, 'seconds...')
        else:
            if (fsize > minsize):
                pass
            else:
                print('FSIZE ERROR JUNK FILE Data not yet available. Waiting', 
                    wait, 'seconds...')
                remove(fpath)
                sleep(wait)
                
                now = datetime.utcnow()
                if ((now-start_time).days >= 1 
                    or (now-start_time).seconds > killtime * 3600):
            
                    exit()

    return params

def gen_paths(init_time):
    from glob import glob
    from os import remove

    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')

    # Purge preexisting index files if any exist and start fresh
    try:
        idxpaths = glob(tmpdir + '%s/models/%s/%s*.idx'%(
            date, ensemble, dateH))
        [remove(idx) for idx in idxpaths]
    except:
        pass
    else:
        print('\nRemoved preexisting index files')

    # Read in the data files [one member at a time]
    member_paths = list()
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

            # Adapt this for init time
            member_paths.append([member, np.sort(glob(tmpdir +
                '%s/models/%s/%s%s*.grib2'%(date, ensemble, dateH, member)))])
        
    return member_paths

def openmfd(paths, lset, cdim):
    """ Dask-free """
    # paths[0] gives the member name
    # paths[1] gives a list of filepaths by hour

    # Open each forecast hour as a dataset...
    datasets = [xr.open_dataset(
            path, engine='cfgrib', #concat_dim=cdim, 
            backend_kwargs={"filter_by_keys":{"typeOfLevel":lset}}) 
                for path in paths[1]]

    if ( (ensemble == 'naefs') & (lset == 'surface') ):
        # Fix the silly issue where CMCE calls tp unknown
        tp = 'unknown' if 'cmce' in str(paths[1]) else 'tp'

        # The NAEFS doesn't include a tp field at hour 0
        # We don't need orog so just swap it for tp with zeros
        datasets[0] = datasets[0].rename({'orog':tp})
        datasets[0][tp].values = np.zeros(datasets[0][tp].shape)

    #...then concatenate them into a single member dataset
    dataset = xr.concat(datasets, dim=cdim)
    dataset = dataset.assign_coords(member_id=paths[0])
    
    if ( ('cmce' in str(paths[1])) & (lset == 'surface') ):
        dataset = dataset.rename({'unknown':'tp'})
    
    # ARW comes in as a total accumulated precip to fhr
    # Deconstruct into 3-hour precip to match the NMB
    if ((lset == 'surface') & ('arw' in paths[0])):
        print("Deconstructing %s from accumulated to step precip"%paths[0])
        arw_tp = np.zeros(dataset.tp.shape)
        arw_tp[0,:,:] = dataset.tp[0,:,:]
        for i in range(1, arw_tp.shape[0]):    
            arw_tp[i,:,:] = dataset['tp'][i,:,:] - dataset['tp'][i-1,:,:]
        
        dataset['tp'].values = arw_tp

    # Clean up the coordinate names a bit
    if lset == 'surface':
        del dataset['surface']
    elif lset == 'isobaricInhPa':
        dataset = dataset.rename({'isobaricInhPa':'level'})
        
    del dataset['step'], dataset['time']

    return dataset

def concat_clean_xarray(xarr_split, cdim):
    
    xarr = xr.concat(xarr_split, dim=cdim)
    xarr = xarr.rename(
        {'number':'member', 'latitude':'lat', 'longitude':'lon'})
    
    xarr['lon'] -= 360

    # Fix member number (otherwise cyclic 0 to mcount)
    xarr.member.values = np.arange(1, xarr.member.size+1)
    
    # Swap the 'valid_time' coordinate for 'time' (which is just init)
    # xarr = xarr.assign_coords(time=xarr.valid_time.values)
    # del xarr['valid_time'], xarr['step']
    xarr = xarr.rename({'valid_time':'time'})
    
    return xarr

def interpolate_prism_daily(doy, year, bounds):
    from netCDF4 import Dataset
    from datetime import datetime
    
    """ Interpolates monthly PRISM totals to a daily total. Assumes the 15th
        (14th for February) is most representative of the month.
        ## Parameters
        doy: The day of year as an int

        year: The year with century as a int

        bounds: A tuple containing boundaries of the domain as indices of
                the PRISM grid. In order xmin, xmax, ymin, ymax.

        prism_dir: The directory the PRISM files are in
        ## Returns
        pclimo: A 2D grid representing a monthly PRISM total were that month
                centered around doy

        Alex Weech
    """

    # Unpack the bounds
    xmin, xmax, ymin, ymax = bounds

    # List of centers of each month
    prism_day = [15] * 12
    prism_day[1] = 14

    # Convert doy and year to a datetime object
    date = datetime.strptime(str(doy) + '-' + str(year), '%j-%Y')

    # Simple case of it being the center day
    center = prism_day[date.month-1]
    if date.day == center:
        prism_path = prism_dir + '/us_' + date.strftime('%m') + '_prcp.nc'
        with Dataset(prism_path, 'r') as prism_cd:
            pclimo = prism_cd.variables['prcp'][0, :, xmin:xmax]
            pclimo = np.flipud(pclimo)[ymin:ymax, :]

    # Else interpolate the two closest months
    else:
        # Check which side of center today is
        if date.day > center:
            month1 = date.month
            year_wrap, month2 = divmod(date.month + 1, 12)
            if month2 == 0:
                year_wrap = 0
                month2 = 12
            centdt1 = datetime(int(year), month1, center)
            centdt2 = datetime(int(year) + year_wrap, month2,
                                  prism_day[month2 - 1])
            weight1 = (date - centdt1).days / (centdt2 - centdt1).days
            weight2 = (centdt2 - date).days / (centdt2 - centdt1).days

        # Else today is before the center
        else:
            month1 = date.month
            year_wrap, month2 = divmod(date.month - 1, 12)
            if month2 == 0:
                year_wrap = -1
                month2 = 12
            centdt1 = datetime(int(year), month1, center)
            centdt2 = datetime(int(year) + year_wrap, month2,
                                  prism_day[month2 - 1])
            weight1 = (centdt1 - date).days / (centdt1 - centdt2).days
            weight2 = (date - centdt2).days / (centdt1 - centdt2).days

        # Open the two files
        file1 = prism_dir + '/us_' + str(month1).zfill(2) + '_prcp.nc'
        file2 = prism_dir + '/us_' + str(month2).zfill(2) + '_prcp.nc'
        with Dataset(file1, 'r') as prism_cd:
            pclimo1 = prism_cd.variables['prcp'][0, :, xmin:xmax]
            pclimo1 = np.flipud(pclimo1)[ymin:ymax, :]

        with Dataset(file2, 'r') as prism_cd:
            pclimo2 = prism_cd.variables['prcp'][0, :, xmin:xmax]
            pclimo2 = np.flipud(pclimo2)[ymin:ymax, :]

        # Interpolate
        pclimo = weight1 * pclimo1 + weight2 * pclimo2

    return pclimo

def downscale_prism(init_time, forecast_time):
    import warnings
    warnings.filterwarnings("ignore")
    
    from scipy import ndimage
    from pandas import to_datetime
    from datetime import datetime, timedelta
    
    # Get the PRISM lats and lons from a sample file
    print('Getting PRISM lats and lons')
    prism = xr.open_dataset(prism_dir + 'us_05_prcp.nc', decode_times=False)

    # Get boundary max and mins using full domain
    xmin = np.max(np.argwhere(prism['lon'].values < -125))
    xmax = np.min(np.argwhere(prism['lon'].values > -100))
    ymin = np.max(np.argwhere(prism['lat'][::-1].values < 30))
    ymax = len(prism['lat'].values) - 1 # Go all the way up
    bounds = (xmin, xmax, ymin, ymax)

    # Subset and mesh
    grid_lons, grid_lats = np.meshgrid(
        prism['lon'][xmin:xmax], prism['lat'][::-1][ymin:ymax])

    # Figure out which days are in this run and put them in a set
    print('Getting PRISM climo')
    date_set = set()

    for i in range(fhrstart, fhrend+1, fhrstep):
        hour_time = init_time + timedelta(hours=i)
        day_of_year = int(hour_time.strftime('%j'))
        date_set.add((day_of_year, hour_time.year))

    # Smoothing algebra
    efold = res_ensemble * 2 / res_prism + 1
    sigma = efold / (np.pi*np.sqrt(2))

    # Loop through the days of this run gathering the climo ratios
    ratios = list()

    for day in date_set:

        pclimo = interpolate_prism_daily(day[0], day[1], bounds)

        # Clip out the missing data
        fixed_prism = np.where(np.logical_and(np.greater(pclimo, 0),
                                              np.isfinite(pclimo)), pclimo, 0)

        # Wyndham's algorithim 
        print('Downscaling PRISM for day of year: {}'.format(
            datetime.strptime(str(day[0]),'%j').strftime('%m/%d')))

        # Create an image smoothed to the model resolution
        smooth_prism = ndimage.filters.gaussian_filter(fixed_prism, sigma,
                                                       mode='nearest')
        smooth_prism = np.where(np.logical_and(np.greater(smooth_prism, 0),
                                               np.isfinite(smooth_prism)),
                                smooth_prism, 0)

        # Divide the real data by the smoothed data to get ratios
        ratios.append([np.where(np.logical_and(np.greater(smooth_prism, 0),
                                               np.greater(fixed_prism, 0)),
                                fixed_prism/smooth_prism, 0), day[0]])
    
    # Sort the prism data back into days (was produced as an unordered set)
    ratios = np.array(ratios)
    ratios = ratios[np.argsort(ratios[:,1].astype(int))]
    
    prism_doy = ratios[:,1]
    prism_data = np.array([x for x in ratios[:,0]])
    
    # Shape into an xarray for easy manipulation
    # Can also save with .to_netcdf if desired
    prism_climo = xr.DataArray(prism_data,
             coords={"time":("time", prism_doy),
                     "lat":(("y", "x"), grid_lats),
                     "lon":(("y", "x"), grid_lons)},
             dims=["time", "y", "x"])
    
    # Do some clipping (Based on Trevor's limits)
    # Not present in old SREF code, added by MW 01/2019
    prism_climo = xr.where(prism_climo < minclip, minclip, prism_climo)
    prism_climo = xr.where(prism_climo > maxclip, maxclip, prism_climo)

    return prism_climo

def get_elev(prism_grid):    
    # Load the elevation DEM
    # Terrainfile is set in config.py
    dem = xr.open_dataset(terrainfile)
    dem = dem.rename({'latitude':'lat', 'longitude':'lon'})

    demlats = dem['lat']
    demlons = dem['lon']
    
    final_lats = prism_grid.lat.values
    final_lons = prism_grid.lon.values

    # As trevor noted, the DEM isn't a perfect match -- 
    # Need to find something better
    xmin = np.where(demlons == demlons.sel(
        lon=final_lons.min(), method='ffill').values)[0][0]
    xmax = np.where(demlons == demlons.sel(
        lon=final_lons.max(), method='bfill').values)[0][0]+1
    ymin = np.where(demlats == demlats.sel(
        lat=final_lats.min(), method='ffill').values)[0][0]
    ymax = np.where(demlats == demlats.sel(
        lat=final_lats.max(), method='bfill').values)[0][0]
    bounds = (xmin, xmax, ymin, ymax)

    elev = dem['elevation'][ymin:ymax, xmin:xmax]
    dem.close()
    
    elevxr = xr.DataArray(elev.values,
         coords={"lat":(("y", "x"), final_lats),
                 "lon":(("y", "x"), final_lons)},
         dims=["y", "x"], name='elev')
    
    return elevxr

def calctw(_t, _rh):
    import warnings
    warnings.filterwarnings("ignore")

    """ Trevor Alcott """
    _tw = (-5.806 + 0.672*(_t-0.006*_t**2 + 
          (0.61 + 0.004*_t + 0.000099*_t**2) * 
          _rh + (-0.000033 - 0.000005*_t - 
          0.0000001*_t**2)*_rh**2))
    
    return xr.where(_tw > _t, _t, _tw)

def calcwbz(_tw, _gh):
    import warnings
    warnings.filterwarnings("ignore")
    from xarray.ufuncs import logical_and
    
    wbz = []
    for i in range(_tw.level.size)[:0:-1]:
        
        # Hi is 'prior' level
        levLO = _tw.level[i-1]
        levHI = _tw.level[i]
        twLO = _tw.isel(level=i-1)
        twHI = _tw.isel(level=i)
        ghLO = _gh.isel(level=i-1)
        ghHI = _gh.isel(level=i)
        
        print('Searching for WBZ between %d and %d hPa'%(levHI, levLO))
        
        twdiff = twLO / (twLO - twHI)
        wbzh = ghLO * twdiff + ghHI * (1 - twdiff)
        
        select = logical_and(twHI < 0., twLO > 0.)
        wbzi = xr.where(select, wbzh, np.nan)
        wbz.append(wbzi)

    return xr.concat(wbz, dim='level').sum(dim='level')

def calct500(_t, _gh, topo):
    
    # Geo Height - Surface Elev + 500 m
    # Gives Geo Heights ABOVE GROUND LEVEL + 500 m buffer
    gh_agl = (_gh - (topo + 500.0)).compute()
    
    # Where this is zero, set to 1.0
    gh_agl = xr.where(gh_agl == 0.0, 1.0, gh_agl)
    
    # If the 1000mb height is > 0, use the 1000 mb temperature to start
    # Otherwise assign t=0
    tvals = xr.where(gh_agl.sel(level=1000) > 0, _t.sel(level=1000), 0) # - 273.15, 0)
    
    for i in range(_t.level.size)[:0:-1]:
        
        # current level
        lc = _t.level.isel(level=i).values
        zc = gh_agl.isel(level=i)
        tc = _t.isel(level=i)# - 273.15
        
        # level above (correct for 'wraparound')
        up = i+1 if i+1 < _t.level.size else 0
        lup = _t.level.isel(level=up).values
        zup = gh_agl.isel(level=up)
        tup = _t.isel(level=up)# - 273.15
        
        # level below
        ldn = _t.level.isel(level=i-1).values
        zdn = gh_agl.isel(level=i-1)
        tdn = _t.isel(level=i-1)# - 273.15
        
        # print(i, lc, lup, ldn)
        
        # Where the geo height AGL > 0 at this level and geo height AGL < 0 at level below...
        tvals = xr.where(((zc > 0.0) & (zdn < 0.0)),
        
        # Do this
        ( ( zc / ( zc - zup ) ) * ( tup - tc ) + tc ),
        
        # Else use tvals already determined
        tvals )
        
    tvals = xr.where(gh_agl.sel(level=500) < 0, _t.sel(level=500), tvals)
        
    return tvals

def calc_slr(t500, wbz, elev):
    ''' Sometimes the old fashioned way of doing things is still the best way.
    Sticking to Trevor's stepwise method which is a little slower but produces
    a reliable result.'''
    
    import warnings
    warnings.filterwarnings("ignore")

    snowlevel = wbz - allsnow
    snowlevel = xr.where(snowlevel < 0., 0., snowlevel)

    initslr = xr.where(t500 < 0., 5. - t500, 5.)
    initslr = xr.where(t500 < -15., 20. + (t500 + 15.), initslr)
    initslr = xr.where(t500 < -20., 15., initslr)

    slr = xr.where(elev >= snowlevel, initslr, 0.)

    slr = xr.where(
        ((elev < snowlevel) & (elev > (snowlevel - melt))),
        (initslr * (elev - (snowlevel - melt)) / melt), slr)

    return slr

def gridfunc(lrdata, lrxy, hrxy):    
    from scipy.interpolate import griddata
    
    hrdata = griddata(lrxy, lrdata.values.flatten(), hrxy, 
                      method='linear', fill_value=0)
    
    hrxr = xr.DataArray(hrdata,
                    coords={"lat":(("y", "x"), hrxy[1]),
                             "lon":(("y", "x"), hrxy[0])},
                    dims=["y", "x"])
    return hrxr

# Do not use... Xarray builds entire netcdf in memory and then dumps
# Extremely memory inefficient.
# def downscale_calc_slr(lr_swapfile, hr_swapfile, iterp_mode='linear'):
#     from os import remove
#     from datetime import datetime
#     from functools import partial
#     from pickle import loads as pl

#     from pandas import to_datetime
    
#     ''' The real meat and potatoes. '''
    
#     lr = pl(np.load(lr_swapfile))
#     hr = pl(np.load(hr_swapfile))

#     dst = datetime.utcnow()
#     mid = lr.member_id.values
    
#     print('Processing member %s'%mid)

#     # Reshape the PRISM array from days to forecast hours
#     forecast_time = lr.time.values
#     forecast_doy = to_datetime(forecast_time).strftime('%j').astype(int)
#     prism_ratios = hr.prism.sel(time=forecast_doy)
#     prism_ratios['time'].values = forecast_time

#     # Set up the low res and hi res xy lat lon arrays
#     if ensemble == 'naefs':
#         lrlon, lrlat = np.meshgrid(lr.lon.values, lr.lat.values)
#         lrxy = (lrlon.flatten(), lrlat.flatten())
#     elif ensemble == 'sref':
#         lrxy = (lr.lon.values.flatten(), lr.lat.values.flatten())
    
#     hrxy = (hr.lon.values, hr.lat.values)
#     gridwrap = partial(gridfunc, lrxy=lrxy, hrxy=hrxy)
    
#     slrst = datetime.utcnow()

#     # Broadcast the elevation to appropriate dimensions...
#     # There's a bug in xr.where() and it fails to do so properly
#     # Submit to github if feeling nice and want to spare others' suffering
#     elev3d = np.broadcast_to(hr.elev.values, (prism_ratios.shape))

#     # Downscale t500, wbz
#     hrt500 = lr.t500.groupby('time').apply(gridwrap)
#     del lr['t500']
#     hrwbz = lr.wbz.groupby('time').apply(gridwrap)
#     del lr['wbz']
#     gc.collect()

#     # print('Downscaled SLR variables for member %s'%mid)
    
#     # Save all vars in dict for easy selection of which to save
#     # Modify in config file only
#     data = {}
#     data['slr'] = calc_slr(hrt500, hrwbz, elev3d)
#     data['slr'] = xr.where(data['slr'] < 0., 0., data['slr'])
#     del hrt500, hrwbz, elev3d
#     gc.collect()
    
#     # print('Calculated SLR for member {} in {}s'.format(
#     #     mid, (datetime.utcnow()-slrst).seconds))
    
#     # Downscale the QPF
#     hrqpf = lr.qpf.groupby('time').apply(gridwrap)
#     data['dqpf'] = hrqpf * prism_ratios.values
#     # print('Downscaled QPF for member %s'%mid)
    
#     del lr['qpf'], prism_ratios, hrqpf, hr
#     gc.collect()
    
#     data['dqpf'] /= 25.4 # mm to inches
#     data['dqpf'] = xr.where(data['dqpf'] < 0., 0., data['dqpf'])
    
#     # Create the hi res downscaled snow grids
#     data['snow'] = data['dqpf'] * data['slr']
#     data['snow'] = xr.where(data['snow'] < 0., 0., data['snow'])
#     # print('Downscaled snow for member %s'%mid)
    
#     # Cumulative sum the dqpf and snow grids to obtain plumes
#     # Easier to save the per-step precip to file and construct later
#     # Rather than save the accumulated. Huge memory drain otherwise.
#     # data['acc_dqpf'] = data['dqpf'].cumsum(dim='time')
#     # data['acc_snow'] = data['snow'].cumsum(dim='time')
    
#     # Set which variables are saved in config.py
#     # print('Saving member %s to netCDF4...'%mid)
#     saveset = xr.Dataset({k:data[k] for k in output_vars})

#     # The metadata is getting lost in the upsample regrid, fix here
#     saveset['member'] = lr.member.values
#     saveset['member_id'] = lr.member_id.values
    
#     # Follow prior directory and naming conventions!
#     inittime = to_datetime(lr.time[0].values)
#     date = inittime.strftime('%Y%m%d')
#     dateH = inittime.strftime('%Y%m%d%H')
    
#     # Write netcdf to temp for speed
#     ncpath = mkdir_p(tmpdir + '%s/models/%s/%s/'%(date, ensemble, dateH))
    
#     filename = '{}_{}_downscaled.nc'.format(dateH, lr.member_id.values)

#     filepath = ncpath + filename
    
#     saveset.to_netcdf(filepath, 
#         format=ncformat)

#     print('Member {} completed at {} in {}s total'.format(
#         mid, datetime.utcnow(), 
#         (datetime.utcnow() - dst).seconds))
    
#     del lr
#     remove(lr_swapfile)

#     return None

def downscale_calc_slr_chunked(lr_swapfile, hr_swapfile, iterp_mode='linear'):
    from os import remove
    from datetime import datetime
    from functools import partial
    from pickle import loads as pl

    from pandas import to_datetime

    ''' The real meat and potatoes. '''
    
    _lr = pl(np.load(lr_swapfile))
    hr = pl(np.load(hr_swapfile))

    dst = datetime.utcnow()
    mid = _lr.member_id.values
    
    inittime = to_datetime(_lr.time[0].values)
    date = inittime.strftime('%Y%m%d')
    dateH = inittime.strftime('%Y%m%d%H')
    tsize = _lr.time.size
    
    print('Processing member %s'%mid)
    
    for i, lr in enumerate(_lr.groupby('time')):
        # print('Processing forecast valid {}'.format(to_datetime(lr[0])))
        forecast_time = [lr[0]]
        lr = lr[1]

        # Reshape the PRISM array from days to forecast hours
        forecast_doy = to_datetime(forecast_time).strftime('%j').astype(int)
        prism_ratios = hr.prism.sel(time=forecast_doy)
        prism_ratios['time'].values = forecast_time

        # Set up the low res and hi res xy lat lon arrays
        if ensemble == 'naefs':
            lrlon, lrlat = np.meshgrid(lr.lon.values, lr.lat.values)
            lrxy = (lrlon.flatten(), lrlat.flatten())
        elif ensemble == 'sref':
            lrxy = (lr.lon.values.flatten(), lr.lat.values.flatten())

        hrxy = (hr.lon.values, hr.lat.values)
        gridwrap = partial(gridfunc, lrxy=lrxy, hrxy=hrxy)
        
        hrt = xr.concat([gridwrap(lr.t.isel(level=l)) 
                         for l in range(lr.level.size)], 
                        dim='level').assign_coords(
                        level = lr.t.level.values)
        
        hrgh = xr.concat([gridwrap(lr.gh.isel(level=l)) 
                         for l in range(lr.level.size)], 
                        dim='level').assign_coords(
                        level = lr.t.level.values)
                
        # Downscale t500, wbz
        hrt500 = calct500(hrt, hrgh, hr.elev)
        del hrt, hrgh, lr['t'], lr['gh'], lr['level']
        gc.collect()
        
        hrwbz = gridwrap(lr.wbz)
        del lr['wbz']
        gc.collect()

        # Save all vars in dict for easy selection of which to save
        # Modify in config file only
        data = {}
        data['slr'] = calc_slr(hrt500, hrwbz, hr.elev.values)
        data['slr'] = xr.where(data['slr'] < 0., 0., data['slr'])
        del hrt500, hrwbz
        gc.collect()

        # Downscale the QPF
        hrqpf = gridwrap(lr.qpf)
        
        prism_ratios = prism_ratios.isel(time=0)
        data['dqpf'] = hrqpf * prism_ratios.values

        del lr['qpf'], prism_ratios, hrqpf
        gc.collect()

        data['dqpf'] /= 25.4 # mm to inches
        data['dqpf'] = xr.where(data['dqpf'] < 0., 0., data['dqpf'])

        # Create the hi res downscaled snow grids
        data['snow'] = data['dqpf'] * data['slr']
        data['snow'] = xr.where(data['snow'] < 0., 0., data['snow'])

        # Cumulative sum the dqpf and snow grids to obtain plumes
        # Easier to save the per-step precip to file and construct later
        # Rather than save the accumulated. Huge memory drain otherwise.
        # data['acc_dqpf'] = data['dqpf'].cumsum(dim='time')
        # data['acc_snow'] = data['snow'].cumsum(dim='time')

        # Set which variables are saved in config.py
        # print('Saving member %s to netCDF4...'%mid)
        saveset = xr.Dataset({k:data[k] for k in output_vars})

        # The metadata is getting lost in the upsample regrid, fix here
        saveset['member'] = lr.member.values
        saveset['member_id'] = lr.member_id.values

        # Follow prior directory and naming conventions!
        # Write netcdf to temp for speed
        ncpath = mkdir_p(tmpdir + '%s/models/%s/%s/'%(date, ensemble, dateH))

        filename = '{}_{}_downscaled.nc'.format(dateH, lr.member_id.values)

        filepath = ncpath + filename
        
        saveset = saveset.expand_dims('time').assign_coords(
            time = forecast_time)
        
        # Write new file or append timestemp depending...
        if to_datetime(forecast_time) == inittime:
            build_netCDF(saveset, i, tsize, dateH, filepath, 'w')
        else:
            build_netCDF(saveset, i, tsize, dateH, filepath, 'a')
        
    print('Member {} completed at {} in {}s total'.format(
        mid, datetime.utcnow(), 
        (datetime.utcnow() - dst).seconds))

    del lr
    remove(lr_swapfile)

    return None

def build_netCDF(xarr, i, tsize, init, fpath, mode):
    from datetime import datetime

    from pandas import to_datetime
    from netCDF4 import Dataset, date2num

    ''' A custom netCDF writer since xarray's doesn't properly append files '''

    with Dataset(fpath, mode, format=ncformat) as dataset:

        if mode == 'w':
            dataset.description = ('Downscaled {} QPF/Snow Grids ' + 
                'Init {} UTC'.format(ensemble.upper(), init))
            dataset.history = 'Created {}'.format(datetime.utcnow())
            dataset.source = 'University of Utah - Steenburgh Research Group'
            
            x = dataset.createDimension('x', xarr.x.size)
            y = dataset.createDimension('y', xarr.y.size)
            t = dataset.createDimension('time', tsize)

            ts = dataset.createVariable('time', 
                'double', ('time',), fill_value=np.nan)
            #ts.calendar = 'gregorian'
            ts.units = 'datetime64' #'hours since 0001-01-01 00:00:00'
            ts.standard_name = 'time'
            ts.long_name = 'time'
            ts.CoordinateAxisType = 'Time'
            
            lat = dataset.createVariable('lat', 
                np.float32, ('y', 'x'), fill_value=np.nan)
            lat.CoordinateAxisType = 'Lat'
            lat.units = 'degrees_north'
            
            lon = dataset.createVariable('lon', 
                np.float32, ('y', 'x'), fill_value=np.nan)
            lon.CoordinateAxisType = 'Lon'
            lon.units = 'degrees_east'
            
            m = dataset.createVariable('member', 'double')
            mid = dataset.createVariable('member_id', 'U10')
                   
            # datenum = date2num(
            #     to_datetime(xarr.time[0].values), 
            #     units=ts.units, calendar=ts.calendar)
            
            ts[i] = xarr.time.values
            lat[:, :] =  xarr.lat.values
            lon[:, :] =  xarr.lon.values
            
            m[0] = xarr.member.values
            mid[:] = xarr.member_id.values

            vardat = {}
            varunits = {'dqpf':'inches', 'snow':'inches', 'slr':'ratio'}
            
            for var in output_vars:
                vardat[var] = dataset.createVariable(
                    var, np.float32, ('time', 'y', 'x'), fill_value=np.nan)
                vardat[var].coordinates = ('lon lat member member_id')
                vardat[var].units = varunits[var]
                
                vardat[var][i, :, :] = xarr[var].values

        elif mode == 'a':
            dataset.variables['time'][i] = xarr.time.values
            for var in output_vars:
                dataset.variables[var][i, : :] = xarr[var].values
            
    return None
        
def dump2swap(xarr, hires, init_time):
    from glob import glob
    from os import path
    from pickle import dumps as pd
    
    print('Dumping to swapfiles on temp as needed')

    # Check if any members already exist, we don't need to waste time 
    # recreating them. 
    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')
    
    tflist = []
    tmppath = mkdir_p(tmpdir + '%s/models/%s/%s/'%(date, ensemble, dateH))
    
    for mno, mid in zip(xarr.member.values, xarr.member_id.values):

        checkfile_archive = tmpdir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
            date, ensemble, dateH, dateH, mid)
        
        # If the file exists in temp, it will be copied later. If the 
        # file already exists in archive, leave it there and the rest will join
        if path.isfile(checkfile_archive):
            pass
        else:
            # Otherwise, save a swapfile for a worker to operate on 
            # and add to list
            tmpfile = tmppath + '%s_swap.npy'%mid
            np.save(tmpfile, pd(xarr.sel(member=mno).compute(), protocol=-1))
            tflist.append(tmpfile)
            
    if len(tflist) > 0:
        hrswapfile = tmppath + 'hires_swap.npy'
        np.save(hrswapfile, pd(hires.compute()))
    else:
        hrswapfile = None
    
    xarr.close()
    del xarr
    gc.collect()
            
    return tflist, hrswapfile

def check_nc_exists(init_time, checkwhere='temp'):
    from os import path, stat
    
    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')
    
    ncfound = []
    
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
            
            if checkwhere == 'temp':
                checkfile = tmpdir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
                    date, ensemble, dateH, dateH, member)
            elif checkwhere == 'archive':
                checkfile = datadir + '%s/models/%s/%s/%s_%s_downscaled.nc'%(
                        date, ensemble, dateH, dateH, member)
            
            # Consider removing this... sometimes an unfinished file 
            # is big enough to pass the test... better to just recreate it
            if path.isfile(checkfile):
                if stat(checkfile).st_size > ncminsize:
                    ncfound.append(True)
                else:
                    ncfound.append(False)
            else:
                ncfound.append(False)
                    
    if (np.where(np.array(ncfound) == False)[0].size) == 0:
        print('Found complete downscaled %s for %s in %s'%(
            ensemble, dateH, checkwhere))

        return True
    
    else:
        return False

def nccopy(ncpath):
    from subprocess import call as sys
    
    cp_command = 'nccopy -d {} {:s} {:s}'.format(
        complevel, ncpath[0], ncpath[1])
    sys(cp_command, shell=True)

    return None
        
def gribcopy(gribpath):
    from subprocess import call as sys
    
    # While I'm aware there is a pythonic copy function, it does
    cp_command = 'cp {:s} {:s}'.format(gribpath[0], gribpath[1])
    sys(cp_command, shell=True)

    return None
    
def temp2archive(init_time, remove_temp=True):
    from glob import glob
    from os import remove as rm
    from multiprocessing import Pool, cpu_count

    dateH = init_time.strftime('%Y%m%d%H')
    date = init_time.strftime('%Y%m%d')

    ncdir = mkdir_p(datadir + '%s/models/%s/%s/'%(
        date, ensemble, dateH))
    nc_temp = glob(tmpdir + '%s/models/%s/%s/%s*.nc'%(
        date, ensemble, dateH, dateH))

    nc_archive = [(ncdir + f.split('/')[-1]) for f in nc_temp]
    nc_paths = [[t, f] for t, f in zip(nc_temp, nc_archive)]
    
    if len(nc_temp) > 0:
        print('Compressing/Copying netCDF files from temp to archive')
        with Pool(cpu_count()-1) as p:
            p.map(nccopy, nc_paths, chunksize=1)
            p.close()
            p.join()
    
    if copy_gribs:
        gribdir = mkdir_p(datadir + '%s/models/%s/'%(date, ensemble))
        grib_temp = glob(tmpdir + '%s/models/%s/%s*.grib2'%(
            date, ensemble, dateH))
        grib_archive = [(gribdir + f.split('/')[-1]) for f in grib_temp]
        grib_paths = [[t, f] for t, f in zip(grib_temp, grib_archive)]
        
        if len(grib_temp) > 0:
            print('Copying grib files from temp to archive')
            with Pool(cpu_count()-1) as p:
                p.map(gribcopy, grib_paths, chunksize=1)
                p.close()
                p.join()

    else:
        grib_temp = glob(tmpdir + '%s/models/%s/%s*.grib2'%(
            date, ensemble, dateH))
        
    # Clean up the temp files if desired
    if remove_temp:
        
        # Note that check_nc_exists returns FALSE if any netcdf files missing
        if check_nc_exists(init_time, checkwhere='archive'):
            print('NetCDF files copied, removing from temp')
            [rm(f) for f in nc_temp]

        else:
            print('Error copying netCDF files, did not delete. Check temp!')

        grib_check = glob(datadir + '%s/models/%s/%s*.grib2'%(
            date, ensemble, dateH))
        
        if len(grib_temp) == len(grib_check):
            print('Removing grib files from temp')
            [rm(f) for f in grib_temp]
        
        grib_idx = glob(tmpdir + '%s/models/%s/%s*.idx'%(
            date, ensemble, dateH))

        [rm(f) for f in grib_idx]
        
    return None