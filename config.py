# Which config file am I?
# Available: ['naefs', 'sref']
ensemble = 'sref'

##########
# mpi
##########

# Number of cores to limit to (overrides all else)
# If None, memory limits still enforced
mpi_limit = None

# Pool worker chunksize
chunksize = 1

# Pool worker spawn method for downscaling jobs
# ['fork', 'forkserver', 'spawn']
spawntype = 'fork'

# Virtual memory needed for pool workers 
# Approx 3.2 GB @ 0-87-3
mem_need = 3.2e9

##########
# paths
##########

# Full chpc path
chpcpath = ('/uufs/chpc.utah.edu/sys/pkg/ldm/' +
            'oper/models/sref/')

# Parent grib/netCDF file directory 
datadir = ('/uufs/chpc.utah.edu/common/home/horel-group/' +
            'archive/')

# Parent grib/netCDF temp file directory
# Set tmpdir = datadir if not using a temp dir
tmpdir = ('/scratch/general/lustre/ldm/rt/')

# Image output directory
imgdir = datadir

# USA terrain file
terrainfile = chpcpath + 'usterrain.nc'

# PRISM files
prism_dir = chpcpath + 'prism/'

# Main script directory
scriptdir = chpcpath + 'scripts/'

# Save the .grib2 files?
copy_gribs = True

###############
# Output files
###############

# Variables to save to netCDF: *(required for plotting)
# {*slr, *dqpf, *snow}
output_vars = ['slr', 'dqpf', 'snow']

# Filetype options
# http://xarray.pydata.org/en/stable/generated/xarray.Dataset.to_netcdf.html
# for notes on compatibility

# Format {‘NETCDF4’, ‘NETCDF4_CLASSIC’, ‘NETCDF3_64BIT’, 'NETCDF3_CLASSIC’}
ncformat = 'NETCDF4'

# Set compression 1-9, with 9 being the smallest file but 
# greatest write time (or use None)
complevel = 9

#################
# ensemble options
#################

# Model init hours to check for
sref_inits = [21, 15, 9, 3]

# Model families
models = ['arw', 'nmb']

# Number of members in each model family
mcount = 13

# Hours to process
# fhrstart must remain 0 but 
# fhrend, fhrstep are mutable
fhrstart = 0
fhrend = 87
fhrstep = 3

# Spatial subset to grab from NOMADS
# use -180 <= lon <= 180
minlon = -130
maxlon = -100
minlat = 30
maxlat = 50

##################
# Download options
##################

# Reject file if it comes in smaller than X bytes (default 10 kb)
minsize = 10000

# When checking for existing ncfile, min size allowed
# 300 mb would still be smaller than using complevel 9
# in almost all cases
ncminsize = 3e8

# Curl command timeout limit (sec)
timeout = 60

# Don't try to download a new run until X
# hours after init time (may use decimal e.g. 4.5)
delay = 4

# wait time between attempts (sec)
#  when data is not yet available
wait = 30

# stop the get function X hours after it starts
# prevents multiple runs downloading at once
killtime = 3

############################
# physics/algorithm options
############################

# Resolution of PRISM in km
res_prism = 0.8

# Resolution of model in km
res_ensemble = 16.0

# Max and min multiplier for QPF downscaling
minclip = 0.3
maxclip = 5.0
    
# SLR Params
allsnow = 150
melt = 300

# # # # # # # # # # # #
# Plotting Options
# # # # # # # # # # # #

# Figures to produce [True or False]
make_maps = True
make_plumes = True

# Scale is relative to Trevor's old figures 
# (made higher res to handle terrain contours and rich colormap)
scale = 1.25

# Plot the elevation contours [True or False]
plot_elevation = False

# (A) How much to smooth the DEM (factor 1-10)
# Low values = smoother terrain, high values = finer detail
# (B) Contour interval, in meters
# Enter dictionary item as {region:(A, B)}
elev_smoothing = {'UT':(3, 1000),
                  'WM':(3, 1000),
                  'CO':(3, 1000),
                  'SN':(3, 1000),
                  'WE':(2, 1000),
                  'NW':(3, 1000),
                  'NU':(6, 500)}

# Thresholds to use for the probability plots (inches)
qpfthresh = [0.01, 1, 2]
snowthresh = [1, 6, 12, 24]

point_locs = {
    'KSLC':(40.77, -111.97), 'CLN':(40.57, -111.65),
    'MTMET':(40.766605, -111.828325), 'PCB':(40.65, -111.52),
    'CLK':(40.67, -111.57), 'PKCU1':(40.6, -111.57),
    'SNI':(41.2, -111.87), 'BLPU1':(41.38, -111.94),
    'TRLU1':(40.68, -110.95), 'DVO':(40.621, -111.496),
    'PWD':(41.37, -111.76), 'SND':(40.37, -111.6),
    'BSK':(45.27, -111.45), 'TOWC2':(40.55, -106.67),
    'CACMP':(40.52, -105.90), 'CABTP':(39.8, -105.77),
    'TUNEP':(39.67, -105.9), 'VAILP':(39.52, -106.22),
    'MONAP':(38.5, -106.32), 'CAMCP':(39.12, -107.27),
    'CARED':(37.9, -107.72), 'CAGMS':(39.05, -108.07),
    'CACBP':(37.7, -107.77), 'CAWCP':(37.47, -106.8),
    'CUMC2':(37.02, -106.45), 'SNO':(47.42, -121.41),
    'CSSL':(39.33, -120.37), 'TGLU1':(41.90, -111.57),
    'SRAC1':(38.32, -119.60), 'SEE':(39.312, -111.429),
    'MSA':(37.65, -119.04), 'MTB':(48.86, -121.68),
    'STV':(47.75, -121.09), 'CMT':(46.93, -121.47),
    'PRS':(46.79, -121.74), 'QUI':(47.47, -123.85),
    'TIM':(45.33, -121.71), 'MBB':(44.0, -121.66),
    'MSB':(41.36, -122.20), 'SQV':(39.2, -120.24),
    'TTSID':(43.86,-114.71), 'JHRV':(43.59, -110.87),
    'LSMU1':(38.48,-109.27),
    #NEW#
    'WDYPK':(40.84148,-111.02062), 'PCCAT':(40.7904494,-111.1053804),
    'STORM':(40.455111, -106.74428), 'SVSID':(43.6612, -114.4121),
    'ASBTP':(35.32581, -111.68559), 'MBNCA':(34.274425, -117.610439),
    'BBEAR':(34.212482, -116.867175), 'TAONM':(36.57443, -105.45328),
    'BASI1':(44.306003, -115.231829), 'LPSI1':(46.63448, -114.58072),
    'BKSI1':(44.62642, -115.79370)
    }

point_names = {
    'KSLC':'Salt Lake City Intl. Airport, UT', 
    'CLN':'Alta Collins, UT',
    'MTMET':'UofU Mountain Met Lab, UT',
    'PCB':'Park City Mountain Resort Base, UT',
    'DVO':'Deer Valley Ontario, UT',
    'SND':'Sundance Mountain Resort, UT',
    'VAILP':'Vail Pass, CO',
    'JHRV':'Jackson Hole Rendezvous Bowl, WY',
    'CSSL':'Central Sierra Snow Lab, CA',
    'SEE':'Seely Creek Manti Skyline, UT',
    'STORM':'Storm Peak Lab',
    'SVSID':'Sun Valley Study Plot',
    'WDYPK':'Windy Peak Uintas',
    'PCCAT':'Park City Powder Cats',
    'PKCU1':'Brighton Crest',
    'BLPU1':'Ben Lomond Peak',
    'ASBTP':'Arizona Snow Bowl Top Patrol',
    'MBNCA':'Mount Baldy Notch CA',
    'BBEAR':'Big Bear CA',
    'TAONM':'Taos Ski Valley, NM',
    'BASI1':'Banner Summit ID',
    'BKSI1':'Big Creek Summit ID',
    'LPSI1':'Lolo Pass ID'
    }

# Region Boundaries (minlat, matlat, minlon, maxlon)
map_regions = {
    #'FULL':(-125, -100, 30, 50),
    'UT':(-114.7, -108.5, 36.7, 42.5),
    'WM':(-117, -108.5, 43, 49),
    'CO':(-110, -104, 36, 41.9),
    'SN':(-123.5, -116.0, 33.5, 41),
    'WE':(-125, -102.5, 31.0, 49.2),
    'NW':(-125, -116.5, 42.0, 49.1),
    'NU':(-113.4, -110.7, 39.5, 41.9),
    }