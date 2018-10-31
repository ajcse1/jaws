from datetime import datetime

import Ngl
import numpy as np
import pandas as pd

import sunposition

try:
    from itertools import izip as zip
except ImportError:  # Python 3.x
    pass


def deg_to_rad(list_deg):
    list_rad = [np.radians(i) for i in list_deg]

    return list_rad


def main(dataset, latitude, longitude, clr_df):
    ddr = 0.25
    rho = 0.8
    smallest_double = 2.2250738585072014e-308

    clrprd_file = clr_df
    clrprd = [(str(x)+'_'+str(y)+'_'+str(z)) for x, y, z in
              zip(clrprd_file[0].tolist(), clrprd_file[1].tolist(), clrprd_file[2].tolist())]

    hours = list(range(24))
    half_hours = list(np.arange(0, 24, 0.5))

    ds = dataset
    df = ds.to_dataframe()

    df.reset_index(level=['time'], inplace=True)
    dates = df['time'].tolist()
    time = sorted(set(dates), key=dates.index)

    tilt_df = pd.DataFrame(index=time, columns=['fsds_adjusted', 'fsds_diff'])

    fsds_adjusted = []
    fsds_diff = []

    lat = latitude
    lon = longitude
    stn = df['station_name'][0]

    grele_path = 'http://grele.ess.uci.edu/jaws/rigb_data/'
    dir_rrtm = 'rrtm-airx3std/'

    for line in clrprd:
        clrdate = line.split('_')[0]
        clrhr_start = int(line.split('_')[1])
        clrhr_end = int(line.split('_')[2])
        year = int(clrdate[:4])
        month = int(clrdate[5:7])
        day = int(clrdate[8:10])
        current_date_hour = []
        for hour in hours:
            current_date_hour.append(datetime(year, month, day, hour))

        calculated_df = pd.DataFrame(index=current_date_hour, columns=['fsds_adjusted', 'fsds_diff'])

        fsds_rrtm = open(grele_path+dir_rrtm+stn+'_'+clrdate.replace('-', '')+'.txt').read().split(',')
        fsds_rrtm = [float(i) for i in fsds_rrtm]

        # Subset dataframe
        df_sub = df[(df.year == year) & (df.month == month) & (df.day == day)]

        sza = df_sub['sza'][0].tolist()
        fsds_jaws = df_sub['shortwave_radiation_down'][0].tolist()
        fsds_jaws_nonmsng = df_sub['shortwave_radiation_down'][0].dropna().tolist()
        indexMissingJAWS = np.where(df_sub['shortwave_radiation_down'][0].isna())
        indexMissingJAWS = [a for b in indexMissingJAWS for a in b]  # Convert to list

        hours_nonmsng = list(range(len(fsds_jaws_nonmsng)))

        # Interpolate fsds and sza for half-hour values
        fsds_intrp = list(Ngl.ftcurv(hours_nonmsng, fsds_jaws_nonmsng, half_hours))

        # Calculate azimuth angle
        az = []
        sza = []
        for hour in hours:
            dtime = datetime(year, month, day, hour, 0)
            az.append(sunposition.sunpos(dtime, lat, lon, 0)[0])
            sza.append(sunposition.sunpos(dtime, lat, lon, 0)[1])
            dtime = datetime(year, month, day, hour, 30)
            az.append(sunposition.sunpos(dtime, lat, lon, 0)[0])
            sza.append(sunposition.sunpos(dtime, lat, lon, 0)[1])

        az = [(i-180) for i in az]

        alpha = [(90-i) for i in sza]

        beta = list(np.arange(0.25, 45.25, 0.25))

        sza_noon = [np.cos(np.radians(i)) for i in sza]

        # Check if measured solar noon time > true solar noon time
        if fsds_intrp.index(max(fsds_intrp)) > sza_noon.index(max(sza_noon)):
            aw = list(np.arange(0, 180, 0.25))
        else:
            aw = list(np.arange(-179.75, 0.25, 0.25))

        az = deg_to_rad(az)
        alpha = deg_to_rad(alpha)
        beta = deg_to_rad(beta)
        aw = deg_to_rad(aw)

        # Make pairs of aw,beta
        pairs = []
        for i in aw:
            for j in beta:
                pairs.append(tuple((i, j)))

        # Find all possible pairs using correct fsds
        possible_pairs = []
        daily_avg_diff = []
        best_pairs = []
        fsds_possiblepair_dict = {}

        for pair in pairs:
            count = 0
            cos_i = []
            fsds_correct = []
            while count < len(alpha):
                cos_i.append((np.cos(alpha[count])*np.cos(az[count]-pair[0])*np.sin(pair[1])+(
                        np.sin(alpha[count])*np.cos(pair[1]))))
                nmr = fsds_intrp[count]*(np.sin(alpha[count])+ddr)
                dnmr = cos_i[count]+(ddr*(1+np.cos(pair[1]))/2.)+(rho*(np.sin(alpha[count])+ddr)*(1-np.cos(pair[1]))/2.)
                if dnmr == 0:
                    dnmr = smallest_double
                fsds_correct.append(nmr/dnmr)
                
                count += 1

            if (abs(cos_i.index(max(cos_i)) - fsds_intrp.index(max(fsds_intrp))) <= 1 and 
               abs(fsds_correct.index(max(fsds_correct)) - sza_noon.index(max(sza_noon))) <= 1):
                possible_pairs.append(pair)

                fsds_correct_half = fsds_correct[::2]
                fsds_possiblepair_dict[pair] = fsds_correct_half

                for msng_idx in indexMissingJAWS:
                    fsds_correct_half.pop(msng_idx)
                    fsds_rrtm.pop(msng_idx)

                diff = [abs(x-y) for x,y in zip(fsds_correct_half[clrhr_start:clrhr_end],
                                                fsds_rrtm[clrhr_start:clrhr_end])]
                daily_avg_diff.append(np.nanmean(diff))


        #########PART-2#############

        dailyavg_possiblepair_dict = dict(zip(daily_avg_diff, possible_pairs))

        if min(dailyavg_possiblepair_dict.keys()) <= 50:
            for val in dailyavg_possiblepair_dict.keys():
                if val <= min(dailyavg_possiblepair_dict.keys())+5:
                    best_pairs.append(dailyavg_possiblepair_dict.get(val))


        #########PART-3#############
 
        fsds_toppair_dict = {k: fsds_possiblepair_dict[k] for k in best_pairs}

        num_spikes = []
        for pair in fsds_toppair_dict:
            fsds_correct_top = fsds_toppair_dict[pair]
            counter = 0
            spike_hrs = 0
            diff_top = [abs(x-y) for x,y in zip(fsds_correct_top, fsds_rrtm)]
            fsds_rrtm_10 = [ij*0.1 for ij in fsds_rrtm]
            for val in diff_top:
                if diff_top[counter] > fsds_rrtm_10[counter]:
                    spike_hrs += 1
                counter += 1
            
            num_spikes.append(spike_hrs)

        top_pair = best_pairs[num_spikes.index(min(num_spikes))]

        fsds_adjusted.append(fsds_toppair_dict[top_pair])
        fsds_diff.append([x-y for x,y in zip(fsds_adjusted, fsds_jaws)])

        calculated_df['fsds_adjusted'] = fsds_adjusted
        calculated_df['fsds_diff'] = fsds_diff

        for val in calculated_df.index:
            tilt_df.at[val, 'fsds_adjusted'] = calculated_df.loc[val]['fsds_adjusted']
            tilt_df.at[val, 'fsds_diff'] = calculated_df.loc[val]['fsds_diff']

    for col in tilt_df:
        tilt_df['fsds_adjusted'] = pd.to_numeric(tilt_df['fsds_adjusted'], errors='coerce')
        tilt_df['fsds_diff'] = pd.to_numeric(tilt_df['fsds_diff'], errors='coerce')

    tilt_df = tilt_df.interpolate()
    fsds_adjusted_values = tilt_df['fsds_adjusted'].tolist()
    fsds_diff_values = tilt_df['fsds_diff'].tolist()

    ds['fsds_adjusted'] = 'time', fsds_adjusted_values
    ds['fsds_diff'] = 'time', fsds_diff_values

    return ds
