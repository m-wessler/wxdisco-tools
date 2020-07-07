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

if __name__ == '__main__':

    site = 'KSLC'

    init_req = sys.argv[1] if len(sys.argv) > 1 else None
    init_time = get_init(init_req)

    pdata = load_plumedata(init_time, site)

    time = pd.to_datetime(pdata.time.values)
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

    timefmt = '%Y-%m-%d %H:%M'
    p1i = np.where(np.array([t.hour for t in time]) == 18)[0][0]

    with open('./test.txt', 'w+') as wfp:
        
        timefmt = '%Y-%m-%d %H:%M'

        wfp.write('{:30s} | {:>30s} | {:>30s} | {:>30s} | {:>30s} |'.format('Init: '+init_time.strftime('%Y-%m-%d %H%M UTC'), 'PoP (%)', 'PQPF [10, 50, 90] (in)', 'PoS (%)', 'PQSF [10, 50, 90] (in)')); wfp.write('\n')
        wfp.write('-'*164); wfp.write('\n')
        wfp.write('{:30s} | '.format('Period 1'))
        wfp.write('{:30s} | '.format('')*4); wfp.write('\n')

        for i in range(11): 

            wfp.write('{:30s} | {:30d} | {:30s} | {:30d} | {:30s} |'.format(
                time[p1i+i].strftime(timefmt),
                pop[p1i+i],
                str(['{:6.3f}'.format(round(v[0], 3)) for v in pqpf.isel(time=[p1i+i]).values]),
                pos[p1i+i],
                str(['{:6.3f}'.format(round(v[0], 3)) for v in pqsf.isel(time=[p1i+i]).values]),
            )); wfp.write('\n')

            if i == 2:
                wfp.write('{:30s} | '.format('')*5); wfp.write('\n')
                wfp.write('{:30s} | '.format('Period 2'))
                wfp.write('{:30s} | '.format('')*4); wfp.write('\n')
                wfp.write('{:30s} | {:30d} | {:30s} | {:30d} | {:30s} |'.format(
                    time[p1i+i].strftime(timefmt),
                    pop[p1i+i],
                    str(['{:6.3f}'.format(round(v[0], 3)) for v in pqpf.isel(time=[p1i+i]).values]),
                    pos[p1i+i],
                    str(['{:6.3f}'.format(round(v[0], 3)) for v in pqsf.isel(time=[p1i+i]).values]),
                )); wfp.write('\n')

    site = 'KSLC'
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
    axs[2].fill_between(time, pqpf.sel(quantile=.1), pqpf.sel(quantile=.9), alpha=0.3, color='g', label='PQPF 5/95 Spread')
    axs[2].set_ylabel('Precipitation [in]\n')
    pqpf_max = np.round(pqpf.sel(quantile=.9).max().values, 2)
    pqpf_max = 1 if pqpf_max < 0.01 else pqpf_max
    axs[2].set_yticks(np.array([np.round(x, 2) for x in np.linspace(0, pqpf_max + (pqpf_max*1.10), 10)]))

    axs[3].plot(time, pqsf.sel(quantile=.5), c='b', label='PQSF Median')
    axs[3].fill_between(time, pqsf.sel(quantile=.1), pqsf.sel(quantile=.9), alpha=0.3, color='b', label='PQSF 5/95 Spread')

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
    plt.savefig('./test.png')
