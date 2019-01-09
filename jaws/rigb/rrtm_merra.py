import csv
from datetime import datetime, timedelta
import os

import climlab
from climlab.domain.field import Field
from climlab.domain import domain
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import xarray as xr

import sunposition

'''
Run RRTM using processed MERRA data from nomiss_merra.py
- run for 8 hour and interpolate for 24 hours
- cal SZA
- output sw_dn
'''


def make_column(lev=None, ps=1013, tmp=None, ts=None):
    state = climlab.column_state(lev=lev)
    num_lev = np.array(lev).size
    lev = state.Tatm.domain.lev

    lev_values = lev.points
    lev_dfs = np.zeros(num_lev)
    lev_dfs[1:] = np.diff(lev_values)
    lev_dfs[0] = lev_values[0]
    lev_dfs = lev_dfs / 2.

    lev.bounds = np.full(num_lev + 1, ps)
    lev.bounds[:-1] = lev_values - lev_dfs
    lev.delta = np.abs(np.diff(lev.bounds))
    sfc, atm = domain.single_column(lev=lev)
    state['Ts'] = Field(ts, domain=sfc)
    state['Tatm'] = Field(tmp, domain=atm)

    return state


def main():
    indir = "nomiss-merra/"
    outdir = "rrtm-merra/"
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
    absorber_vmr['CFC11'] = 0.280 / 1E9
    absorber_vmr['CFC12'] = 0.503 / 1E9
    absorber_vmr['CFC22'] = 0.
    absorber_vmr['CCL4'] = 0.

    # stn names and lat/lon
    lst_stn = pd.read_csv('../resources/stations_radiation.txt')
    stn_names = lst_stn['network_name'].tolist()
    latstn = lst_stn['lat'].tolist()
    lonstn = lst_stn['lon'].tolist()

    nstn = len(stn_names)
    '''stn_names = 'gcnet_summit'
    latstn = 72.57972
    lonstn = -38.50454
    nstn = 1'''

    start_year = 2002
    end_year = 2003
    years = [str(i) for i in list(range(start_year, end_year, 1))]
    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    days = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
            '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']

    hours_merra = [1, 4, 7, 10, 13, 16, 19, 22]
    hours_24 = list(range(24))

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
                    sw_dn = []

                    flname = indir + stn + '.' + str(year) + str(month) + str(day) + '.nc'

                    if not os.path.isfile(flname):
                        continue

                    print(flname)
                    fl_date = flname.split('.')[1]

                    fin = xr.open_dataset(flname)
                    fout = outdir + stn + '.' + fl_date + '.txt'

                    tmp = fin['t'].values
                    ts = fin['ts'].values
                    plev = fin['plev'].values
                    ps = fin['ps'].values
                    h2o_q = fin['q'].values
                    o3_mmr = fin['o3'].values
                    o3_vmr = climlab.utils.thermo.mmr_to_vmr(o3_mmr, gas='O3')

                    aod_count = fin['aod_count'].values

                    # knob
                    aod = np.zeros((6, 1, 72))
                    aod[1, 0, -aod_count:] = 0.12 / aod_count
                    aod[5, 0, :15] = 0.0077 / 15

                    idx = 0

                    for hr in hours_merra:
                        dtime = datetime.strptime(fl_date, "%Y%m%d") + timedelta(hours=hr, minutes=30)

                        sza = sunposition.sunpos(dtime, lat_deg, lon_deg, 0)[1]
                        cossza = np.cos(np.radians(sza))

                        state = make_column(lev=plev[idx], ps=ps, tmp=tmp[idx], ts=ts)
                        absorber_vmr['O3'] = o3_vmr[idx]

                        rad = climlab.radiation.RRTMG(name='Radiation', state=state, specific_humidity=h2o_q[idx],
                                                      albedo=alb, coszen=cossza, absorber_vmr=absorber_vmr,
                                                      emissivity=emis, S0=solar_constant, icld=0, iaer=6,
                                                      ecaer_sw=aod)
                        rad.compute_diagnostics()

                        dout = rad.to_xarray(diagnostics=True)
                        sw_dn.append(dout['SW_flux_down_clr'].values[-1])

                        idx += 1

                    sw_dn_24 = CubicSpline(hours_merra, sw_dn, extrapolate=True)(hours_24)

                    with open(fout, 'w') as outfile:
                        wr = csv.writer(outfile)
                        wr.writerow(sw_dn_24)


if __name__ == '__main__':
    main()
