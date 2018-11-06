import csv
from datetime import datetime, timedelta
import glob
import os

import climlab
from climlab.domain.field import Field
from climlab.domain import domain
import numpy as np
import pandas as pd
import xarray as xr

import sunposition

'''
Run RRTM using processed AIRS data from nomiss_airs.py
- run 24 hour twice (A and D)
- cal SZA
- output netCDF
'''


def make_column(lev=None,  ps=1013, tmp=None, ts=None):
    state = climlab.column_state(lev=lev)
    num_lev = np.array(lev).size
    lev = state.Tatm.domain.lev

    lev_values = lev.points
    lev_dfs = np.zeros(num_lev)
    lev_dfs[1:] = np.diff(lev_values)
    lev_dfs[0] = lev_values[0]
    lev_dfs = lev_dfs/2.

    lev.bounds = np.full(num_lev+1, ps)
    lev.bounds[:-1] = lev_values-lev_dfs
    lev.delta = np.abs(np.diff(lev.bounds))
    sfc, atm = domain.single_column(lev=lev)
    state['Ts'] = Field(ts, domain=sfc)
    state['Tatm'] = Field(tmp, domain=atm)
    
    return state


def main():
    indir = "nomiss-airx3std/"
    outdir = "rrtm-airx3std/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # constants
    alb = 0.75
    emis = 0.985
    solar_constant = 1367.0

    absorber_vmr = dict()
    absorber_vmr['CO2'] = 355. / 1E6
    absorber_vmr['CH4'] = 1714. / 1E9
    absorber_vmr['N2O'] = 311. / 1E9
    absorber_vmr['O2'] = 0.21
    absorber_vmr['CFC11'] = 0.280/1E9
    absorber_vmr['CFC12'] = 0.503/1E9
    absorber_vmr['CFC22'] = 0.
    absorber_vmr['CCL4'] = 0.

    # stn names and lat/lon
    lst_stn = pd.read_csv('stations_radiation.txt')
    stn_names = lst_stn['stn_name'].tolist()
    latstn = lst_stn['lat'].tolist()
    lonstn = lst_stn['lon'].tolist()

    nstn = len(stn_names)
    '''stn_names = 'imau09'
    latstn = -75.0
    lonstn = 0.0
    nstn = 1'''

    start_year = 2010
    end_year = 2017
    years = [str(i) for i in list(range(start_year, end_year, 1))]
    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
            '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']

    # main function
    for i in range(nstn):
        stn = stn_names[i]
        lat_deg = latstn[i]
        lon_deg = lonstn[i]
        '''stn = stn_names
        lat_deg = latstn
        lon_deg = lonstn'''

        for year in years:
            for month in months:
                for day in days:
                    try:
                        sw_dn = []
                        sw_dn_final = [None]*24
                        for sfx in ['A', 'D']:
                            fn = indir + stn + '_' + str(year) + str(month) + str(day) + '_' + sfx + '.nc'

                            fin = xr.open_dataset(fn)

                            tmp = fin['t'].values
                            # ts = fin['sfc_tmp'].values
                            ts = fin['ts'].values
                            plev = fin['plev'].values
                            # ps = fin['sfc_prs'].values
                            ps = fin['ps'].values

                            state = make_column(lev=plev, ps=ps, tmp=tmp, ts=ts)

                            o3 = fin['o3'].values
                            absorber_vmr['O3'] = o3

                            h2o_q = fin['q'].values

                            aod_count = fin['aod_count'].values

                            # knob
                            aod = np.zeros((6, 1, 24))
                            aod[1, 0, -aod_count:] = 0.12 / aod_count
                            aod[5, 0, :15] = 0.0077 / 15

                            for hr in range(24):
                                dtime = datetime.strptime(fn.split('_')[1], "%Y%m%d") + timedelta(hours=hr)

                                # fo = outdir+stn+'_'+dtime.strftime('%Y%m%d:%H')+fn[-5:]
                                fo = outdir+stn+'_'+dtime.strftime('%Y%m%d')+'.txt'

                                sza = sunposition.sunpos(dtime, lat_deg, lon_deg, 0)[1]
                                cossza = np.cos(np.radians(sza))

                                rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o_q,
                                                              albedo=alb, coszen=cossza, absorber_vmr=absorber_vmr,
                                                              emissivity=emis, S0=solar_constant, icld=0, aod=aod)
                                rad.compute_diagnostics()

                                dout = rad.to_xarray(diagnostics=True)
                                sw_dn.append(dout['SW_flux_down_clr'].values[-1])
                                # break
                                # dout.to_netcdf(fo)

                        count = 0
                        while count < 24:
                            # sw_dn_final[count] = [(g+h)/2.0 for g,h in zip(sw_dn[count], sw_dn[count+24])]
                            sw_dn_final[count] = (sw_dn[count]+sw_dn[count+24])/2.0
                            count += 1

                        with open(fo, 'w') as outfile:
                            wr = csv.writer(outfile)
                            wr.writerow(sw_dn_final)
                        # print(fo)
                        # np.savetxt(fo, sw_dn_final, delimiter=",", fmt='%s')

                        # outfile = open(fo, 'w')
                        # myString = ",".join(map(str,sw_dn_final))
                        # outfile.write("%s," % myString)

                        # for item in sw_dn_final:
                        #    outfile.write("%s," % item)

                        # break
                    except:
                        pass


if __name__ == '__main__':
    main()
