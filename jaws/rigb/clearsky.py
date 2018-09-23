#import warnings
#warnings.filterwarnings('ignore')
import pandas as pd
import xarray as xr
import numpy as np
import scipy.interpolate

from itertools import groupby
from operator import itemgetter

from datetime import datetime, timedelta

import Ngl

import matplotlib.pyplot as plt


def interpolate(x,y):
    dv = scipy.interpolate.interp1d(x,y, fill_value='extrapolate')
    alpha = 0.01
    dv_left, dv_right, tg_left, tg_right = ([None]*len(x) for i in range(4))
    i = 0
    for hour in x:
        dv_left[i] = float(dv(i-alpha))
        dv_right[i] = float(dv(i+alpha))
        tg_left[i] = (y[i] - dv_left[i])/alpha
        tg_right[i] = (dv_right[i] - y[i])/alpha
        i += 1
    diff = [x1 - x2 for (x1,x2) in zip(tg_left, tg_right)]
    return diff


def clr_prd(dat_sza, tg_fsds, tg_sza, hours, out_file, year, day_of_year, count):
    scale = 1400
    offset = 0
    offset_range = 15

    tg_sza_scale = [i*scale for i in tg_sza]
    tg_sza_up = [i-offset+offset_range for i in tg_sza_scale]
    tg_sza_dn = [i-offset-offset_range for i in tg_sza_scale]

    clr_hrs = []
    daylight = []
    for hour in hours:
        if dat_sza[hour] > 0:
            daylight.append(hour)

        if (hour==0):
            clr_hrs.append(hour)
        elif (hour>0) and (hour<23):
            if (tg_fsds[hour] < tg_sza_up[hour]) and (tg_fsds[hour] > tg_sza_dn[hour]):
                clr_hrs.append(hour)
        elif (hour==23):
            clr_hrs.append(hour)
    
    cons_clr_hrs = []
    for k, g in groupby(enumerate(clr_hrs), lambda ix : ix[0] - ix[1]):
        cons_clr_hrs.append(list(map(itemgetter(1), g)))
    
    final_hrs = write_to_file(cons_clr_hrs, daylight, out_file, year, day_of_year, count)
    
    return final_hrs


def clr_shift(final_hrs, dat_sza, dat_fill, hrs):
    if not final_hrs:
        dat_sza_shift = list(range(len(dat_sza)))

        shift = [1,2,3]

        for i in shift:
            if dat_sza.index(max(dat_sza)) < dat_fill.argmax:
                dat_sza_shift[-i:] = [0] * (len(dat_sza_shift) - i)
                dat_sza_shift[:-i] = dat_sza[i:]

            elif dat_sza.index(max(dat_sza)) > dat_fill.argmax:
                dat_sza_shift[:i] = [0] * (len(dat_sza_shift) - i)
                dat_sza_shift[i:] = dat_sza[:-i]

            tg_sza_shift = interpolate(hrs,dat_sza_shift)

            cons_clr_hrs = clr_prd(dat_sza, tg_fsds, tg_sza_shift, hours, out_file, year, day_of_year, count)

            if cons_clr_hrs:
                break


def write_to_file(cons_clr_hrs, daylight, out_file, year, day_of_year, count):
    final_hrs = []
    for group in cons_clr_hrs:
        if len(group) >= int(len(daylight)/2):
            final_hrs.append(group)
            
            with open(out_file, "a") as fl:
                fl.write("{}, {}, {}\n".format(
                    datetime.strptime("%s %s" % (year,int(day_of_year[count])), "%Y %j").date(), group[0], group[-1]))
            return final_hrs


#def main(args):
def main():
    global tg_fsds, dat_sza, dat_fill, hrs, hours, out_file, year, day_of_year, count

    #infile = args.input_file
    infile = 'UPE-L.nc'
    out_file = "cleardays.txt"

    with open(out_file, "w") as fl:
        pass


    ds = xr.open_dataset(infile)
    ds = ds.drop('time_bounds')
    df = ds.to_dataframe()


    year = 2010
    df_sub = df[(df.year==year) & (df.month>4) & (df.month<9)]
    #df_sub = df[(df.year==year) & (df.month==7)]


    day_of_year = df_sub['day_of_year'].tolist()
    day_of_year = list(set(day_of_year))

    months = list(range(5,9))
    days = list(range(1,32))
    hours = list(range(24))

    #plt.figure(figsize=(20,40))
    #plot_number = 1
    count = 0

    for month in months:
        for day in days:
            try:
                df_temp = df_sub[(df_sub.month==month) & (df_sub.day==day)]
                
                dat = df_temp['shortwave_radiation_down'].tolist()
                dat_rmvmsng = df_temp['shortwave_radiation_down'].dropna().tolist()
                hrs = list(range(len(dat)))
                hrs_rmvmsng = list(range(len(dat_rmvmsng)))
                dat_fill = Ngl.ftcurv(hrs_rmvmsng, dat_rmvmsng, hrs)
                dat_sza = [np.cos(np.radians(i)) for i in df_temp['sza'].tolist()]

                tg_fsds = interpolate(hrs,dat_fill)
                tg_sza = interpolate(hrs,dat_sza)
                
                final_hrs = clr_prd(dat_sza, tg_fsds, tg_sza, hours, out_file, year, day_of_year, count)
                
                cons_clr_hrs = clr_shift(final_hrs, dat_sza, dat_fill, hrs)
                
                
                count += 1

                '''ax1 = plt.subplot(31, 4, plot_number)
                ax1.plot(hours, dat_fill)

                ax2 = ax1.twinx()
                ax2.plot(hours, dat_sza, linestyle='--', color='hotpink')

                if day in [4,8,12,16,20,24,28,31]:
                    pass
                else:
                    ax2.set_yticks([])

                plot_number += 1'''
            except:
                continue


#if __name__ == '__main__':
#    main()