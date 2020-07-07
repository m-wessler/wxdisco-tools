#!/uufs/chpc.utah.edu/sys/installdir/anaconda3/2018.12/bin/python

import sys
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from datetime import datetime

from funcs import *
from config import *
from plotfuncs import *

import warnings
warnings.filterwarnings('ignore')

wd = '/uufs/chpc.utah.edu/common/home/u1070830/public_html/wxdisco/'

if __name__ == '__main__':

    site = 'KSLC'

    init_req = sys.argv[1] if len(sys.argv) > 1 else None
    init_time = get_init(init_req)

    pdata = load_plumedata(init_time, site)

    time = pd.to_datetime(pdata.time.values)
    timefmt = '%Y-%m-%d %H:%M'
    p1i = np.where(np.array([t.hour for t in time]) == 18)[0][0]

    pop, pos = [], []
    for i, t in enumerate(time):
        n = pdata.isel(time=i).dqpf.size

        p_yes = np.where(pdata.isel(time=i).dqpf >= 0.01)[0].size
        pop.append(np.round(p_yes/n*100, 0).astype(int))

        s_yes = np.where(pdata.isel(time=i).snow >= 0.1)[0].size
        pos.append(np.round(s_yes/n*100, 0).astype(int))
        
    pop, pos = np.array(pop), np.array(pos)
    pqpf = pdata.dqpf.quantile([.1, .5, .9], dim='member')
    pqsf = pdata.snow.quantile([.1, .5, .9], dim='member')

    p1 = pdata.isel(time=slice(p1i+1, p1i+3))
    n1 = p1.isel(time=0).member.size

    p1_pop = int(np.round((np.where(p1.dqpf.sum(dim='time') > 0.01)[0].size/n1)*100, 0))
    p1_pos = int(np.round((np.where(p1.snow.sum(dim='time') > 0.1)[0].size/n1)*100, 0))

    p2 = pdata.isel(time=slice(p1i+4, p1i+11))
    n2 = p2.isel(time=0).member.size

    p2_pop = int(np.round((np.where(p2.dqpf.sum(dim='time') > 0.01)[0].size/n2)*100, 0))
    p2_pos = int(np.round((np.where(p2.snow.sum(dim='time') > 0.1)[0].size/n2)*100, 0))

    p1['dqpf'].values[np.where(p1.dqpf <= 0.01)] = 0.
    p1['snow'].values[np.where(p1.dqpf <= 0.1)] = 0.
    p1_pqpf = [np.round(v, 2) for v in p1.dqpf.sum(dim='time').quantile([.1, .5, .9], dim='member').values]
    p1_pqsf = [np.round(v, 2) for v in p1.snow.sum(dim='time').quantile([.1, .5, .9], dim='member').values]

    p2['dqpf'].values[np.where(p2.dqpf <= 0.01)] = 0.
    p2['snow'].values[np.where(p2.dqpf <= 0.1)] = 0.
    p2_pqpf = p2.dqpf.sum(dim='time').quantile([.1, .5, .9], dim='member').values
    p2_pqsf = p2.snow.sum(dim='time').quantile([.1, .5, .9], dim='member').values

    try:
        with open(wd + 'forecast.csv', 'r+') as rfp:
            lines = rfp.readlines()
            lenlines = len(lines)
    except:
        lines = []
        lenlines = 0

    maxlines = 91 if lenlines > 91 else None

    with open(wd + 'forecast.csv', 'w+') as wfp:
        wfp.write(' , , , , \n')
        
        wfp.write('{}, {}, {}, {}, {}\n'.format(
            'KSLC Downscaled SREF Init: '+init_time.strftime('%Y-%m-%d %H%M UTC'), 
            'PoP (%)', 
            'PQPF 10/50/90 (in)', 
            'PoS (%)', 
            'PQSF 10/50/90 (in)'))

        for i in range(11): 
            # wfp.write('{}, {}, {}, {}, {}\n'.format(
            #     '3-h period ending: ' + time[p1i+i].strftime(timefmt),
            #     pop[p1i+i],
            #     str(['{:6.3f}'.format(round(v[0], 3)) for v in pqpf.isel(time=[p1i+i]).values]).replace(',', ' ').replace("'",'').replace('[', '').replace(']', ''),
            #     pos[p1i+i],
            #     str(['{:6.3f}'.format(round(v[0], 3)) for v in pqsf.isel(time=[p1i+i]).values]).replace(',', ' ').replace("'",'').replace('[', '').replace(']', '')))
            
            if i in [2, 10]:
                period = 'Period 1 (%s – %s)'%(time[p1i], time[p1i+i]) if i == 2 else 'Period 2 (%s – %s)'%(time[p1i+10], time[p1i+i])
                _pop, _pos = (p1_pop, p1_pos) if i == 2 else (p2_pop, p2_pos)
                _pqpf, _pqsf = (p1_pqpf, p1_pqsf) if i == 2 else (p2_pqpf, p2_pqsf)
                
                wfp.write('{}, {}, {}, {}, {}\n'.format(
                    period,
                    _pop,
                    str(['{:6.3f}'.format(round(v, 3)) for v in _pqpf]).replace(',', ' ').replace("'",'').replace('[', '').replace(']', ''),
                    _pos,
                    str(['{:6.3f}'.format(round(v, 3)) for v in _pqsf]).replace(',', ' ').replace("'",'').replace('[', '').replace(']', '')
                ))

        [wfp.write(line) for line in lines[:maxlines]]

    plt.rcParams.update({'font.size': 16})

    fig, axs = plt.subplots(2, 2, figsize=(30, 16), facecolor='w', sharex=True)
    fig.subplots_adjust(hspace=0.03, wspace=0.015)
    axs = axs.flatten()

    axs[0].scatter(time, pop, s=250, marker="+", color='g', label='PoP')
    axs[0].set_ylabel('Probability of Precipitation\n')
    axs[0].set_yticks(np.arange(0, 101, 10))

    axs[1].scatter(time, pos, s=250, marker="+", color='b', label='PoS')
    axs[1].set_ylim(top=100)
    axs[1].set_yticks(np.arange(0, 101, 10))

    axs[2].plot(time, pqpf.sel(quantile=.5), c='g', label='PQPF Median')
    axs[2].fill_between(time, pqpf.sel(quantile=.1), pqpf.sel(quantile=.9), alpha=0.3, color='g', label='PQPF 10/90 Spread')
    axs[2].set_ylabel('Precipitation [in]\n')
    pqpf_max = np.round(pqpf.sel(quantile=.9).max().values, 2)
    pqpf_max = 1 if pqpf_max < 0.01 else pqpf_max
    axs[2].set_yticks(np.array([np.round(x, 2) for x in np.linspace(0, pqpf_max + (pqpf_max*1.10), 10)]))

    axs[3].plot(time, pqsf.sel(quantile=.5), c='b', label='PQSF Median')
    axs[3].fill_between(time, pqsf.sel(quantile=.1), pqsf.sel(quantile=.9), alpha=0.3, color='b', label='PQSF 10/90 Spread')

    pqsf_max = np.round(pqsf.sel(quantile=.9).max().values, 2)
    pqsf_max = 1 if pqsf_max < 0.1 else pqsf_max
    axs[3].set_yticks(np.array([np.round(x, 2) for x in np.linspace(0, pqsf_max + (pqsf_max*1.10), 10)]))

    for i, ax in enumerate(axs):    
        if i in [1, 3]:
            ax.yaxis.tick_right()
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper right', fontsize='x-large')
        ax.set_xlim([time[0], time[-1]])
        ax.set_xlabel('\nTime [UTC]')

        hours = mdates.HourLocator(interval = 6)
        h_fmt = mdates.DateFormatter('%m/%d %H:%M')
        ax.xaxis.set_major_locator(hours)
        ax.xaxis.set_major_formatter(h_fmt)
        
        for tick in ax.get_xticklabels():
                tick.set_rotation(70)
            
        if i in [0, 1]:
            ax.set_xticklabels([])   

        ax.grid(True)
            
    plt.suptitle('%s Downscaled SREF Precipitation Probabiltiies\nInitialized %s UTC\n' %(site, init_time), y=.95, x=0.51)
    plt.savefig(wd + 'current_sref.png')
