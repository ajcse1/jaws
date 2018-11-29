from datetime import datetime
import os
import requests

import numpy as np
import pandas as pd
import xarray as xr

try:
    from jaws import common
except:
    import common


def deg_to_rad(list_deg):
    list_rad = [np.radians(i) for i in list_deg]

    return list_rad


def main(dataset, args):
    rho = 0.8  # reflectance
    smallest_double = 2.2250738585072014e-308
    fillvalue_double = 9.969209968386869e+36
    dtime_1970, tz = common.time_common(args.tz)
    idx_count = 0

    ds = dataset
    ds = ds.drop('time_bounds')
    df = ds.to_dataframe()

    date_hour = [datetime.fromtimestamp(i, tz) for i in df.index.values]
    dates = [i.date() for i in date_hour]
    dates = sorted(set(dates), key=dates.index)

    df['fsds_corrected'] = ''

    df.reset_index(level=['time'], inplace=True)
    stn_name = df['station_name'][0]

    grele_path = 'http://grele.ess.uci.edu/jaws/rigb_data/'
    dir_ceres = 'ceres/'
    sfx = '.ceres.nc'
    url = grele_path + dir_ceres + stn_name + sfx
    r = requests.get(url, allow_redirects=True)
    open(stn_name + sfx, 'wb').write(r.content)

    for date in dates:
        year = date.year
        month = date.month
        day = date.day

        ceres_df = xr.open_dataset(stn_name + sfx).to_dataframe()
        cf = ceres_df.loc[str(year)+'-'+str(month)+'-'+str(day):str(year)+'-'+str(month)+'-'+str(day)]['cldarea_total_1h'].values.tolist()
        cf = [i/100 for i in cf]

        df_sub = df[(df.year == year) & (df.month == month) & (df.day == day)]

        fsds_jaws = df_sub['sw_down'].tolist()
        fsds_jaws = [fillvalue_double if np.isnan(i) else i for i in fsds_jaws]

        sza = df_sub['sza'].tolist()
        az = df_sub['az'].tolist()

        az = [(i - 180) for i in az]

        alpha = [(90 - i) for i in sza]

        aw = df_sub['tilt_direction'].tolist()

        beta = df_sub['tilt_angle'].tolist()

        az = deg_to_rad(az)
        alpha = deg_to_rad(alpha)
        beta = deg_to_rad(beta)
        aw = deg_to_rad(aw)

        count = 0
        while count < len(alpha):
            cos_i = (np.cos(alpha[count]) * np.cos(az[count] - aw[count]) * np.sin(beta[count]) + (
                    np.sin(alpha[count]) * np.cos(beta[count])))

            ddr = (0.2+0.8*cf[count])/(0.8-0.8*cf[count])

            nmr = fsds_jaws[count] * (np.sin(alpha[count]) + ddr)

            dnmr = cos_i + (ddr * (1 + np.cos(beta[count])) / 2.) + (
                        rho * (np.sin(alpha[count]) + ddr) * (1 - np.cos(beta[count])) / 2.)
            if dnmr == 0:
                dnmr = smallest_double

            df.at[idx_count, 'fsds_corrected'] = nmr/dnmr

            count += 1
            idx_count += 1

    fsds_corrected_values = df['fsds_corrected'].tolist()
    dataset['fsds_corrected'] = 'time', fsds_corrected_values

    os.remove(stn_name + sfx)

    return dataset
