import obspy
from obspy import UTCDateTime
from obspy.core.event import *
from obspy.core.event.header import *
import pandas as pd
import numpy as np

from geopy.distance import great_circle

'''ANSSDF2CAT converts ANSS df to ObsPy Catalog'''
def anssdf2cat(df):
	cat = Catalog()
	for i, row in df.iterrows():
    
        #print(i)
        #print(row)
        
		new_origin = Origin()
		new_origin.time = UTCDateTime(row['DateTime'])
		new_origin.latitude = float(row['Latitude'])
		new_origin.longitude = float(row['Longitude'])
		new_origin.depth = float(row['Depth'])
        
        # only populate Magnitude object if magnitude is reported
		new_magnitude = Magnitude()
		if ~(np.isnan(float(row['Magnitude']))):
			new_magnitude.mag = float(row['Magnitude'])
			new_magnitude.magnitude_type = row['MagType']
        #new_magnitude.evaluation_mode = EvaluationMode('automatic')
        
		new_event = Event()
		new_event.origins = [new_origin]
		new_event.magnitudes = [new_magnitude]
        
		cat.append(new_event)
        #print('')

	print(' ')
	print(cat)
	return cat

'''ANSSCSV2CAT converts an ANSS csv file to ObsPy Catalog'''
def ansscsv2cat(file):
	import pandas as pd
	df = pd.read_csv(file)
	cat = anssdf2cat(df)
	return cat


'''RFILTER Radial filter of lat/lon values based on an annulus in meters'''
'''
Examples:
>>> rfilter(vlatlon, eqlatlon, r)
where vlatlon is either a NumPy array or a tuple
where r is either a NumPy array (of size 1 or 2) or a tuple or an Int
where eqlatlon is always a n-by-2 NumPy array

>>> eqlatlon1 = rfilter(vlatlon, eqlatlon, r)
'''
def rfilter(vlatlon, eqlatlon, r):
	# convert radial input into tuple
	if (type(r)==type(np.array([]))):
		if r.size==2:
			print('Annulus filter is a NumPy array of size 2')
			r = tuple(r)
		elif r.size==1:
			print('Annulus filter is NumPy array of size 1')
			r = (0, r[0])
	elif type(r)==type((0,0)):
		#print('Annulus filter is a Tuple')
		r = r
	elif type(r)==type(0):
		#print('Annulus filter is an Int')
		r = (0,r)
	else:
		print('Radius must be a tuple or an Integer')

	r_km = np.array(list(map(lambda i: great_circle(vlatlon, i).meters, eqlatlon))) # radial distance from volcano to eq
	return eqlatlon[(r_km > r[0]) & (r_km <= r[1])]