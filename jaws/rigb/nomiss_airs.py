import xarray as xr
import numpy as np
import pandas as pd
import Ngl
import glob
import os
import warnings
warnings.filterwarnings('ignore')


def getDSLatLon(infile):
    ds = xr.open_dataset(infile)

    vars = ['StdPressureLev:ascending_TqJoint', 'SurfPres_Forecast_TqJ_A', 'SurfPres_Forecast_TqJ_D', 
            'Temperature_TqJ_A', 'Temperature_TqJ_D', 'SurfAirTemp_TqJ_A', 'SurfAirTemp_TqJ_D',
            'SurfSkinTemp_TqJ_A', 'SurfSkinTemp_TqJ_D', 'H2O_MMR_TqJ_A', 'H2O_MMR_TqJ_D', 
            'H2O_MMR_Surf_TqJ_A', 'H2O_MMR_Surf_TqJ_D', 'O3_VMR_TqJ_A', 'O3_VMR_TqJ_D']
    ds_sub = ds[vars]

    var_lat = ds_sub['YDim:ascending_TqJoint'].values
    var_lon = ds_sub['XDim:ascending_TqJoint'].values

    lat_lon_var = []

    for i in var_lat:
        for j in var_lon:
            lat_lon_var.append([i,j])

    return ds_sub, lat_lon_var


def haversine_np(latlon1, latlon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    """
    latlon1, latlon2 = map(np.radians, [latlon1, latlon2])
    
    i = 0
    stn_new=[]
    
    for item in latlon1:
        j=0
        dist=[]
        for value in latlon2:
            dlat = latlon1[i][0]-latlon2[j][0]
            dlon = latlon1[i][1]-latlon2[j][1]
            
            a = np.sin(dlat/2.0)**2 + np.cos(latlon2[j][0]) * np.cos(latlon1[i][0]) * np.sin(dlon/2.0)**2
            c = 2 * np.arcsin(np.sqrt(a))
            km = 6367 * c
            dist.append(km)
            j += 1
        idx=dist.index(min(dist))
        stn_new.append(latlon2[idx])
        i += 1
    
    stn_new = list(map(np.degrees, stn_new))
    return stn_new


def getStnLatLonName(station_file):
    lst_stn = pd.read_csv(station_file)
    stn_names = lst_stn['stn_name'].tolist()

    latstn = lst_stn['lat'].tolist()
    lonstn = lst_stn['lon'].tolist()

    '''stn_names = 'imau09'
    latstn = -75.0
    lonstn = 0.0'''

    lat_lon_stn = np.column_stack((latstn, lonstn))

    return stn_names, lat_lon_stn

def main():
    #indir = './'
    indir = '../../../../data/wenshanw/airs/hdf_airx3std/'
    outdir = 'nomiss-airx3std'
    nplev = 24
    nplev_h2o = 12
    pin, tin, qin = ([None]*(nplev+1) for i in range(3))
    sfxs = ['A', 'D']
    q_top = 1e-7
    stn_names, lat_lon_stn = getStnLatLonName('../resources/stations_radiation.txt')

    #for infile in glob.iglob(indir+'*.hdf'):
    for infile in glob.iglob(indir+'AIRS.2009.12.*.hdf'):
        print(infile)
        i=0
        ds_sub, lat_lon_var = getDSLatLon(infile)
        stn_new = haversine_np(lat_lon_stn,lat_lon_var)
        
        for item in stn_new:
            try:
                j = 0
                
                temp = {'YDim:ascending_TqJoint': round(stn_new[i][j],1), 'YDim:descending_TqJoint': round(stn_new[i][j],1),
                        'XDim:ascending_TqJoint': round(stn_new[i][j+1],1), 'XDim:descending_TqJoint': round(stn_new[i][j+1],1)
                       }
                ds_sub_temp = ds_sub.sel(temp)
                
                ####################################
                plev = ds_sub_temp['StdPressureLev:ascending_TqJoint'].values.tolist()[::-1]

                for sfx in sfxs:
                    #Pressure
                    pin[:nplev] = plev[:]
                    v_ps_airs = 'SurfPres_Forecast_TqJ_'+sfx
                    ps_airs = ds_sub_temp[v_ps_airs].values.tolist()
                    if ps_airs:
                        pin[-1] = ps_airs

                    #Temperature
                    v_t_airs = 'Temperature_TqJ_'+sfx
                    tin[:nplev] = ds_sub_temp[v_t_airs].values.tolist()[::-1]
                    v_ta_airs = 'SurfAirTemp_TqJ_'+sfx
                    ta_airs = ds_sub_temp[v_ta_airs].values.tolist()
                    if ta_airs:
                        tin[-1] = ta_airs

                    v_ts_airs = 'SurfSkinTemp_TqJ_'+sfx
                    ts_airs = ds_sub_temp[v_ts_airs].values.tolist()
                    #Water vapor mixing ratio
                    v_q_airs = 'H2O_MMR_TqJ_'+sfx
                    qin[0] = q_top
                    qin[(nplev-nplev_h2o):nplev] = [alpha/1000. for alpha in ds_sub_temp[v_q_airs].values.tolist()[::-1]]

                    v_qs_airs = 'H2O_MMR_Surf_TqJ_'+sfx
                    qs_airs = ds_sub_temp[v_qs_airs].values.tolist()
                    if qs_airs:
                        qin[-1] = qs_airs/1000.

                    # All pressure >= ps
                    pin_sub = pin[:-1]

                    x=0
                    idx_lst=[]
                    for val in pin_sub:
                            if val > (pin[-1]-5):
                                idx_lst.append(pin_sub.index(pin_sub[x]))
                            x += 1
                    
                    for idx in idx_lst:
                        for lst in [pin,tin,qin]:
                            lst[idx] = None
                    
                    #Fill missing
                    pin = pd.Series(pin)
                    pin = pin.interpolate(limit_direction='both').values.tolist()
                    pout = pin[:-1]

                    def log_interpolation(var_name, **kwargs):
                        var_in = np.array(var_name, dtype=np.float64)
                        nans, x = np.isnan(var_in), lambda z: z.nonzero()[0]
                        if var_name == tin:
                            var_in[nans] = np.exp(np.interp(np.log(x(nans)), np.log(x(~nans)), np.log(var_in[~nans])))

                            var_out = var_in[:-1]
                        
                        elif var_name == qin:
                            var_in[nans] = np.interp(x(nans), x(~nans), var_in[~nans])
                            var_out = var_in[:-1]
                        elif var_name == o3_airs:
                            '''pout_rmvmsng = pout
                            nan_idx = np.argwhere(np.isnan(np.array(o3_airs)))
                            nan_idx = [ijk for ijk in nan_idx]
                            for idx in sorted(nan_idx, reverse=True):
                                del pout_rmvmsng[int(idx)]
                    
                            print(len(pout_rmvmsng))
                            print(len(var_in[~nans]))'''
                            
                            var_in[nans] = np.interp(x(nans), x(~nans), var_in[~nans])
                            #var_in = np.exp(np.interp(np.log(plev), np.log(pout_rmvmsng), np.log(var_in[~nans])))
                            #var_in = np.exp(np.interp(np.log(plev), np.log(pout), np.log(var_in)))
                            var_out = var_in
                        
                        return var_out

                    try:
                        #print qin
                        tout = log_interpolation(tin)
                        qout = log_interpolation(qin)

                        #Convert q_mmr to q_specific-humidity
                        qout_q = [beta/(beta+1) for beta in qout]
                    except:
                        tout = tin[:-1]
                        qout_q = qin[:-1]

                    #Ozone
                    v_o3_airs = 'O3_VMR_TqJ_'+sfx
                    o3_airs = ds_sub_temp[v_o3_airs].values.tolist()[::-1]

                    try:
                        o3_out = log_interpolation(o3_airs, plev = plev, pout = pout)
                    except:
                        o3_out = o3_airs
                    #print(o3_out)


                    ds = xr.Dataset({'plev': ('PLEV', pout), 't': ('PLEV', tout), 'q': ('PLEV', qout_q), 'o3': ('PLEV', o3_out), 
                            'ts': ts_airs, 'ta': tin[-1], 'ps': ps_airs})
                    
                    ####################################
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)
                    basename = os.path.basename(infile)
                    outfile = outdir+'/'+stn_names[i]+'_'+basename[5:15].replace('.','')+'_'+sfx+'.nc'

                    ds.to_netcdf(outfile)
            
            except:
                pass
            
            #break

            i += 1


if __name__ == '__main__':
    main()