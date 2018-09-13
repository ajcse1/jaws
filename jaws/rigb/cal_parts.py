import numpy as np
import pandas as pd
import xarray as xr
import Ngl
import sunposition
from datetime import datetime
from math import sin, cos
import csv


def main():
	ds = xr.open_dataset('promice_Upernavik-L_20090817_20170916.nc')
	df = ds.to_dataframe()

	latitude = df['latitude'][0][0]
	longitude = df['longitude'][0][0]

	ddr = 0.25
	pho = 0.8

	year = 2010
	month = 7
	days = [2,3,5,6,7,8,9,10,12,17]
	hours = range(24)
	half_hours = list(np.arange(0,24,0.5))


	for day in days:
		
		#Subset dataframe
		df_sub = df[(df.year == year) & (df.month == month) & (df.day == day)]

		fsds = df_sub['shortwave_radiation_down'][0].tolist()
		sza = df_sub['sza'][0].tolist()

		#Interpolate fsds and sza for half-hour values
		fsds_intrp = list(Ngl.ftcurv(hours, fsds, half_hours))
		sza_intrp = list(Ngl.ftcurv(hours, sza, half_hours))

		
		#Calculate azimuth angle
		az=[]
		for hour in hours:
		    dtime = datetime(year,month,day,hour,0)
		    az.append(sunposition.sunpos(dtime, latitude, longitude, 0)[0])

		az_intrp = list(Ngl.ftcurv(hours, az, half_hours))

		alpha = [90-i for i in sza_intrp]

		beta = list(np.arange(0,45,0.25))

		
		#Check if measured solar noon time > true solar noon time
		if fsds_intrp.index(max(fsds_intrp)) > sza_intrp.index(max(sza_intrp)):
		    aw = list(np.arange(0,180,0.25))
		else:
		    aw = list(np.arange(-179.75,0.25,0.25))


		#Make pairs of aw,beta
		pairs = []
		for i in aw:
		    for j in beta:
		        pairs.append(tuple((i,j)))


		#Find all possible pairs using correct fsds
		possible_pairs = []

		for pair in pairs:
		    count = 0
		    cos_i = []
		    correct_fsds = []
		    while count < len(alpha):
		        cos_i.append((cos(alpha[count])*cos(az_intrp[count]-pair[0])*sin(pair[1])+(sin(alpha[count])*cos(pair[1]))))
		        nmr = fsds_intrp[count]*(sin(alpha[count])+ddr)
		        dnmr = cos_i[count]+(ddr*(1+cos(pair[1]))/2.)+(pho*(sin(alpha[count])+ddr)*(1-cos(pair[1]))/2.)
		        correct_fsds.append(nmr/dnmr)
		        
		        count += 1

		    if (abs(cos_i.index(max(cos_i)) - fsds_intrp.index(max(fsds_intrp))) <= 0.5 and 
		       abs(correct_fsds.index(max(correct_fsds)) - sza_intrp.index(max(sza_intrp))) <= 0.5):
		        possible_pairs.append(pair)


		print(day)
		print(possible_pairs)
		print('***********')

if __name__ == '__main__':
	main()