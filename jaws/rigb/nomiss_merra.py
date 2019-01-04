import xarray as xr
import numpy as np
import pandas as pd
import glob
import os
import warnings

warnings.filterwarnings('ignore')


def get_ds_latlon(infile):
    ds = xr.open_dataset(infile)

    var_lat = ds['lat'].values
    var_lon = ds['lon'].values

    lat_lon_var = []

    for i in var_lat:
        for j in var_lon:
            lat_lon_var.append([i, j])

    return ds, lat_lon_var


def haversine_np(latlon1, latlon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    """
    latlon1, latlon2 = map(np.radians, [latlon1, latlon2])

    i = 0
    stn_new = []

    for item in latlon1:
        j = 0
        dist = []
        for value in latlon2:
            dlat = latlon1[i][0] - latlon2[j][0]
            dlon = latlon1[i][1] - latlon2[j][1]

            a = np.sin(dlat / 2.0) ** 2 + np.cos(latlon2[j][0]) * np.cos(latlon1[i][0]) * np.sin(dlon / 2.0) ** 2
            c = 2 * np.arcsin(np.sqrt(a))
            km = 6367 * c
            dist.append(km)
            j += 1
        idx = dist.index(min(dist))
        stn_new.append(latlon2[idx])
        i += 1

    stn_new = list(map(np.degrees, stn_new))
    return stn_new


def get_stn_latlonname(station_file):
    lst_stn = pd.read_csv(station_file)
    stn_names = lst_stn['network_name'].tolist()

    latstn = lst_stn['lat'].tolist()
    lonstn = lst_stn['lon'].tolist()
    lonstn = [360 + i if i < 0 else i for i in lonstn]

    lat_lon_stn = np.column_stack((latstn, lonstn))

    return stn_names, lat_lon_stn


def log_interpolation(var_name, **kwargs):
    var_in = np.array(var_name, dtype=np.float64)
    nans, x = np.isnan(var_in), lambda z: z.nonzero()[0]
    if var_name == tin:
        var_in[nans] = np.exp(np.interp(np.log(x(nans)), np.log(x(~nans)), np.log(var_in[~nans])))

        var_out = var_in[:-1]

    elif var_name == qin:
        var_in[nans] = np.interp(x(nans), x(~nans), var_in[~nans])
        var_out = var_in[:-1]
    elif var_name == o3_merra:
        '''pout_rmvmsng = pout
        nan_idx = np.argwhere(np.isnan(np.array(o3_merra)))
        nan_idx = [ijk for ijk in nan_idx]
        for idx in sorted(nan_idx, reverse=True):
            del pout_rmvmsng[int(idx)]

        print(len(pout_rmvmsng))
        print(len(var_in[~nans]))'''

        var_in[nans] = np.interp(x(nans), x(~nans), var_in[~nans])
        # var_in = np.exp(np.interp(np.log(plev), np.log(pout_rmvmsng), np.log(var_in[~nans])))
        # var_in = np.exp(np.interp(np.log(plev), np.log(pout), np.log(var_in)))
        var_out = var_in

    return var_out


def checkna(variable):
    variable = pd.to_numeric(variable, errors='coerce')
    if np.count_nonzero(~np.isnan(variable)) >= len(variable) / 3:
        return True
    else:
        return False


def main():
    # indir = '/data/wenshanw/merra2/rigb_input/'
    indir = 'C:/Users/Ajay/data/merra2/'
    outdir = 'merra_nomiss-airx3std'

    nplev = 72

    global tin, qin, o3_merra
    pin, tin, qin = ([None] * (nplev + 1) for _ in range(3))

    stn_names, lat_lon_stn = get_stn_latlonname('../resources/stations_radiation.txt')

    for infile in glob.iglob(indir + '*20020701*.nc'):

        print(infile)

        x_coord = 0
        ds, lat_lon_var = get_ds_latlon(infile)
        stn_new = haversine_np(lat_lon_stn, lat_lon_var)

        for item in stn_new:
            y_coord = 0

            temp = {'lat': round(stn_new[x_coord][y_coord], 1),
                    'lon': round(stn_new[x_coord][y_coord + 1], 1)
                    }
            ds_temp = ds.sel(temp)  # Subset the dataset for only these dimensions

            # Pressure
            plev = ds_temp['PL'].values.tolist()
            plev = [i*100 for i in plev]
            pin[:nplev] = plev[:]

            ps_merra = ds_temp['PS'].values.tolist()
            ps_merra = [i*100 for i in ps_merra]
            if not np.isnan(ps_merra):
                pin[-1] = ps_merra

            # Temperature
            tin[:nplev] = ds_temp['T'].values.tolist()
            ta_merra = ds_temp['TLML'].values.tolist()
            if not np.isnan(ta_merra):
                tin[-1] = ta_merra

            ts_merra = ta_merra

            # Water vapor mixing ratio
            qin[:nplev] = ds_temp['QV'].values.tolist()
            qs_merra = ds_temp['QLML'].values.tolist()
            if not np.isnan(qs_merra):
                qin[-1] = qs_merra

            # Ozone
            o3_merra = ds_temp['O3'].values.tolist()

            # Get index of all values, where pressure >= surface pressure
            pressureSafetyMeasure = 5
            pin_sub = pin[:-1]  # All values except surface pressure
            x = 0
            idx_lst = []
            for val in pin_sub:
                if val > (pin[-1] - pressureSafetyMeasure):
                    idx_lst.append(pin_sub.index(pin_sub[x]))
                x += 1

            # Set the values at above indexes as missing
            for idx in idx_lst:
                for lst in [pin, tin, qin]:
                    lst[idx] = None

            # Interpolate the above missing values
            pin = pd.Series(pin)
            pin = pin.interpolate(limit_direction='both').values.tolist()
            pout = pin[:-1]

            aod_count = 0
            for press in pout:
                if (pin[-1] - 100) < press < pin[-1]:
                    aod_count += 1

            # If more than (2/3)rd values are missing for any variable, skip that day
            if checkna(tin) and checkna(qin) and checkna(o3_merra):
                tout = log_interpolation(tin)
                qout = log_interpolation(qin)
                o3_out = log_interpolation(o3_merra, plev=plev, pout=pout)
            else:
                x_coord += 1

            # Write variables to xarray-dataset
            ds = xr.Dataset({'plev': ('PLEV', pout), 't': ('PLEV', tout), 'q': ('PLEV', qout),
                             'o3': ('PLEV', o3_out), 'ts': ts_merra, 'ta': tin[-1], 'ps': ps_merra,
                             'aod_count': aod_count})

            # Write to netCDF file
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            basename = os.path.basename(infile)
            outfile = outdir + '/' + stn_names[x_coord] + '.' + basename[5:15].replace('.', '') + '.nc'

            ds.to_netcdf(outfile)

            x_coord += 1


if __name__ == '__main__':
    main()
