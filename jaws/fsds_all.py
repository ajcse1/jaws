from datetime import datetime, timedelta
import requests

import Ngl
import numpy as np
import pandas as pd
import xarray as xr

import sunposition

try:
    from itertools import izip as zip
except ImportError:  # Python 3.x
    pass


def deg_to_rad(list_deg):
    list_rad = [np.radians(i) for i in list_deg]

    return list_rad


def main(dataset):
    rho = 0.8  # reflectance
    smallest_double = 2.2250738585072014e-308
    fillvalue_double = 9.969209968386869e+36

    ds = dataset
    ds = ds.drop('time_bounds')
    df = ds.to_dataframe()

    df['corrected_fsds'] = ''

    df.reset_index(level=['time'], inplace=True)
    dates = df['time'].tolist()
    #time = sorted(set(dates), key=dates.index)

    stn = df['station_name'][0]
    year = df['year'].tolist()
    day_of_year = df['day_of_year'].tolist()

    grele_path = 'http://grele.ess.uci.edu/jaws/rigb_data/'
    dir_ceres = 'ceres/'
    url = grele_path + dir_ceres + '.ceres.nc'
    r = requests.get(url, allow_redirects=True)
    open(stn + '.ceres.nc', 'wb').write(r.content)

    for yr in year:
        for doy in day_of_year:
            dt = (datetime(yr, 1, 1) + timedelta(days=doy))
            month = dt.month
            day = dt.day
            day -= 1  # Day of year starts with 1 instead of 0

            cf_ds = xr.open_dataset(stn + '.' + yr + '.ceres.nc')
            cf_df = cf_ds.to_dataframe()
            cf = cf_df.loc[yr+'-'+month+'-'+day:yr+'-'+month+'-'+day]['cldarea_total_1h'].values.tolist()
            cf = [i/100 for i in cf]

            df_sub = df[(df.year == year) & (df.month == month) & (df.day == day)]

            fsds_jaws = df_sub['shortwave_radiation_down'][0].tolist()
            fsds_jaws = [fillvalue_double if i.isna() else i for i in fsds_jaws]

            sza = df_sub['sza'][0].tolist()
            az = df_sub['az'][0].tolist()

            az = [(i - 180) for i in az]

            alpha = [(90 - i) for i in sza]

            aw = df_sub['tilt_direction'][0].tolist()

            beta = df_sub['tilt_angle'][0].tolist()

            az = deg_to_rad(az)
            alpha = deg_to_rad(alpha)
            beta = deg_to_rad(beta)
            aw = deg_to_rad(aw)

            count = 0
            fsds_correct = []
            while count < len(alpha):
                cos_i = (np.cos(alpha[count]) * np.cos(az[count] - aw) * np.sin(beta) + (
                        np.sin(alpha[count]) * np.cos(beta)))

                ddr = (0.2+0.8*cf[count])/(0.8-0.8*cf[count])

                nmr = fsds_jaws[count] * (np.sin(alpha[count]) + ddr)

                dnmr = cos_i + (ddr * (1 + np.cos(beta)) / 2.) + (
                            rho * (np.sin(alpha[count]) + ddr) * (1 - np.cos(beta)) / 2.)
                if dnmr == 0:
                    dnmr = smallest_double

                fsds_correct.append(nmr / dnmr)
                df.at[count_outer, 'fsds_corrected'] = nmr/dnmr

                count += 1
