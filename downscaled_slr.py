import gc
import os
import psutil
import numpy as np
from sys import argv
from functools import partial
from datetime import datetime
from multiprocessing import Pool, cpu_count, get_context

import xarray as xr

from funcs import *
from config import *
from plotfuncs import *

import warnings
warnings.filterwarnings('ignore')

os.environ['OMP_NUM_THREADS'] = '1'

if __name__ == '__main__':
    scriptstart = datetime.utcnow()

    # Acquire input time if any
    init_req = argv[1] if len(argv) == 2 else None

    # Otherwise determine init time
    init_time = get_init(init_req)
    print('Retrieving %s forecast initalized %s'%(ensemble, init_time))
    
    # Make sure we don't needlessly produce files that already exist
    # Helpful for repeat plotting
    temp_exists = check_nc_exists(init_time, checkwhere='temp')
    archive_exists = check_nc_exists(init_time, checkwhere='archive')

    if archive_exists:
        usetemp = False
        process_new = False

    elif temp_exists:
        usetemp = True
        process_new = False

    else:
        usetemp = True
        process_new = True

    if process_new:
        # Do the downloading
        get_grids(init_time)

        # Gather paths to member files
        mpaths = gen_paths(init_time)

        # Read members without dask distributed (slow)
        print('\nReading individual member files')

        cores = len(mpaths) if len(mpaths) <= cpu_count() else cpu_count()
        mpi_limit = mpi_limit if mpi_limit is not None else 999
        cores = cores if cores <= mpi_limit else mpi_limit

        with get_context(spawntype).Pool(cores) as p:
            openmfd_wrap = partial(openmfd, lset='surface', cdim='valid_time')
            _surface = p.map(openmfd_wrap, mpaths)

            openmfd_wrap = partial(
                openmfd, lset='isobaricInhPa', cdim='valid_time')
            _isobaric = p.map(openmfd_wrap, mpaths)
            
            p.close()
            p.join()

        # Gather to single xarray
        print('Concatenating members...', end='')
        surface = concat_clean_xarray(_surface, cdim='number')
        isobaric = concat_clean_xarray(_isobaric, cdim='number')

        [f.close() for f in _surface]
        [f.close() for f in _isobaric]
        del _surface, _isobaric
        print('...done')

        # Downscale prism to establish the ratios and output grid
        print('\nGetting PRISM ratios...')
        prism = downscale_prism(init_time, surface.time.values)

        # Get elevation data from the DEM and clip to the PRISM grid
        print('\nReading elevation data...')
        elev = get_elev(prism)

        # Start with calculating tw at each isobaric level
        gh = isobaric['gh']
        t = isobaric['t'] - 273.15
        rh = isobaric['r']

        # Calculate tw(3D), wbz(2D)
        print('\nCalculating wet bulb temp... ', end='')
        tw = calctw(t, rh)

        print(' ...done\nCalculating height of wet bulb zero...')
        wbz = calcwbz(tw, gh)

        # Free up memory and nuke no longer needed variables
        del rh, tw
        gc.collect()
        print(' ...done\n')

        # Low Res
        qpf = surface.tp

        del isobaric, surface
        gc.collect()

        # Since float32 is sufficient, reduce type to save on memory
        lowres = xr.Dataset({'qpf':qpf.astype(np.float32), 
                            't':t.astype(np.float32),
                            'gh':gh.astype(np.float32),
                            'wbz':wbz.astype(np.float32)})
        lowres['lat'] = lowres['lat'].astype(np.float32)
        lowres['lon'] = lowres['lon'].astype(np.float32)

        # Get rid of unnecessary tagalongs
        del qpf

        # Repack the high-resolution data into a neat dataset
        hires = xr.Dataset({'prism':prism, 'elev':elev})

        # Break the members into a list to iterate over for Python Pool Workers
        lowres_swapfiles, hires_swapfile = dump2swap(lowres, hires, init_time)

        # Close grib files, clean up memory
        lowres.close()
        hires.close()
        del lowres, hires
        gc.collect()

        # Here we are! Nice and speedy
        print('Low Res Data Ready (Completed in {})\n'.format(
            datetime.utcnow() - scriptstart))

        dstart = datetime.now()

        # Wrap the function and arguments with partial
        downscale_mpi = partial(
                            downscale_calc_slr_chunked, 
                            hr_swapfile=hires_swapfile, 
                            iterp_mode='linear')

        # Figure out how many pool workers we can get away with
        mem_avail = psutil.virtual_memory().available
        # Using ceil, if causing memory problems in the future, use floor
        memlim = np.floor(mem_avail/mem_need).astype(int)

        # Friendly, readable printouts
        hmem_avail = bytes2human(mem_avail) 
        hmem_need = bytes2human(mem_need)

        # Set the limits - CORES is memory limited automatically!
        cores = (memlim if memlim <= len(lowres_swapfiles) 
            else len(lowres_swapfiles))
        mpi_limit = mpi_limit if mpi_limit is not None else 9999
        cores = cores if cores <= mpi_limit else mpi_limit
        
        # Hard cutoff - allows the SREF in one go but limits the NAEFS
        # Due to diminishing returns with increase cpu usage
        # cores = cores if cores <= 26 else 26

        print(("\nMem Avail: {}\nMem Needed: {}\nWorkers Needed: {}\n" +
                "Cores Available: {}\nCores to use: {}\n").format(
                hmem_avail, hmem_need, len(lowres_swapfiles), 
                memlim, cores))

        print('Starting %d pool workers'%cores)
        # Using a container ensures nice and neat garbage collection
        with get_context(spawntype).Pool(cores) as p:
            p.map(downscale_mpi, lowres_swapfiles, chunksize=chunksize)
            p.close()
            p.join()

        os.remove(hires_swapfile)

        print('\nData processing completed at {} in: {}'.format(
            datetime.now(), datetime.utcnow() - scriptstart))
    
    if make_maps == True:
        # Make the plots - First the spatial data (large in memory, 
        # so load once and distribute only critical arrays as dict)
        print('Loading QPF data for spatial plots')
        qpf_acc = load_mapdata(init_time, 'dqpf', temp=usetemp)
        print('Calculating QPF statistics for spatial plots')
        qpf_mapdata = calc_stats(qpf_acc, init_time)
        del qpf_acc; gc.collect()

        print('Loading Snow data for spatial plots')
        snow_acc = load_mapdata(init_time, 'snow', temp=usetemp)
        print('Calculating Snow statistics for spatial plots')
        snow_mapdata = calc_stats(snow_acc, init_time)
        snow_mapdata['mqpf'] = qpf_mapdata['mean']
        del snow_acc; gc.collect()

        qpfplot = partial(ensemble_qpf6panel, dd=qpf_mapdata)
        snowplot = partial(ensemble_snow6panel, dd=snow_mapdata)

        nmaps = len(map_regions.keys())*2
        cpus = cpu_count() - 1
        cores = nmaps if nmaps <= cpus else cpus
        
        # Forcing 'spawn' for plots due to data mismatch issues with fork
        with get_context('spawn').Pool(cores) as p:
            print('Producing %d spatial plots on %d cores'%(
                nmaps, cores))

            p.map(qpfplot, map_regions.keys(), chunksize=1)
            p.map(snowplot, map_regions.keys(), chunksize=1)

            p.close()
            p.join()

        # Free up some memory
        del qpf_mapdata, snow_mapdata
        gc.collect()

    if make_plumes == True:
        plumeplot = partial(ensemble_plume, init_time=init_time, temp=usetemp)
        plumelocs = point_locs.keys()

        nplumes = len(plumelocs)
        cpus = cpu_count() - 1
        cores = nplumes if nplumes <= cpus else cpus

        # Forcing 'spawn' for plots due to data mismatch issues with fork
        with get_context('spawn').Pool(cores) as p:
            print('Producing %d plumes on %d cores'%(
                len(plumelocs), cores))

            p.map(plumeplot, plumelocs, chunksize=1)
            p.close()
            p.join()

    print('\nDownscaling and plotting completed at {} UTC in: {}'.format(
            datetime.utcnow(), datetime.utcnow() - scriptstart))

    if usetemp == True:
        temp2archive(init_time)

    print('\nScript Complete at {} UTC in {}\n'.format(
    datetime.utcnow(), datetime.utcnow() - scriptstart))

    # Despite calling close() in functions some h5netcdf files are
    # lingering open... force them to close before exiting
    # for obj in gc.get_objects():
    #     try:
    #         obj.close()
    #     except:
    #         pass

