# Airwave alarm to be run on set of channels in small aperture infrasound array
# Based on MATLAB code originally written by Matt Haney and John Lyons
#
# Wech 2017-06-08

from obspy import UTCDateTime, Stream
from obspy.core.util import AttribDict
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.cross_correlation import xcorr
import numpy as np
from os import remove
from pandas import DataFrame
import time


### This code probably doesn't work yet
def add_coordinate_info(st,SCNL):
    '''Adds scnl info to traces
    Parameters
    ----------
    SCNL=[
        {'scnl':'ADKI.BDF.AV.01'	, 'sta_lat': 51.86190727	, 'sta_lon': -176.6438598},
        {'scnl':'ADKI.BDF.AV.02'	, 'sta_lat': 51.86324162	, 'sta_lon': -176.6435998},
        {'scnl':'ADKI.BDF.AV.03'	, 'sta_lat': 51.86226962	, 'sta_lon': -176.6446503},
        {'scnl':'ADKI.BDF.AV.04'	, 'sta_lat': 51.86246609	, 'sta_lon': -176.6457851},
        {'scnl':'ADKI.BDF.AV.05'	, 'sta_lat': 51.86326916	, 'sta_lon': -176.6461231},
        {'scnl':'ADKI.BDF.AV.06'	, 'sta_lat': 51.86157572	, 'sta_lon': -176.6469340},
    ]
    '''
    #### compare remaining stations with lat/lon station info in config file
    #### to attach lat/lon info with each corresponding trace


    for tr in st:
        if tr.stats.location=='':
            tr.stats.location='--'
        tmp_scnl='{}.{}.{}.{}'.format(tr.stats.station,
                                      tr.stats.channel,
                                      tr.stats.network,
                                      tr.stats.location)
        tmp_lat=SCNL[SCNL['scnl']==tmp_scnl].sta_lat.values[0]
        tmp_lon=SCNL[SCNL['scnl']==tmp_scnl].sta_lon.values[0]
        tr.stats.coordinates=AttribDict({
                                'latitude': tmp_lat,
                                'longitude': tmp_lon,
                                'elevation': 0.0})
    return st


def get_volcano_backazimuth(st, vlatlon):
    '''Returns the backazimuth
    
    Parameters
    -----------
    st          : ObsPy Stream
        each trace must contain latlon info, e.g.,
            tr.stats.coordinates.longitude and
            tr.stats.coordinates.latitude
    latlon      : tuple
        latitude and longitude of target
        
    Returns
    -------
    az[1]      :
        backazimuth
    '''
    lon0=np.mean([tr.stats.coordinates.longitude for tr in st])
    lat0=np.mean([tr.stats.coordinates.latitude for tr in st])
    az=gps2dist_azimuth(lat0,lon0,vlatlon[0],vlatlon[1])
    return az[1]


def setup_coordinate_system(st):
    '''Returns a coordinate system for traces with latlon info
    Each trace must contain tr.stats.coordinates.longitude, tr.stats.coordinates.latitude
    '''
    R = 6372.7976   # radius of the earth
    lons  = np.array([tr.stats.coordinates.longitude for tr in st])
    lats  = np.array([tr.stats.coordinates.latitude for tr in st])
    lon0  = lons.mean()*np.pi/180.0
    lat0  = lats.mean()*np.pi/180.0
    yx    = R*np.array([ lats*np.pi/180.0-lat0, (lons*np.pi/180.0-lon0)*np.cos(lat0) ]).T
    intsd = np.zeros([len(lons),len(lons)])
    ints_az= np.zeros([len(lons),len(lons)])
    for ii in range(len(st[:-1])):
        for jj in range(ii+1,len(st)):
            # intsd[i,j]=np.sqrt(np.square(yx[j][0]-yx[i][0])+np.square(yx[j][1]-yx[i][1]))
            tmp=gps2dist_azimuth(lats[ii],lons[ii],lats[jj],lons[jj])
            intsd[ii,jj]=tmp[0]
            ints_az[ii,jj]=tmp[1]

    return yx, intsd, ints_az


def calc_triggers(st, intsd, cc_shift_length=3*50, min_cc=0.6, vmin=0.28):
    '''
    cc_shift_length     : float
        maximum samples to shift in cross-correlation (usually at 50 sps)
    min_cc              : float
        min normalized correlation coefficient to accept
    vmin                : float
        minimum velocity to accept
    
    '''
    lags       = np.array([])
    lags_inds1 = np.array([])
    lags_inds2 = np.array([])
    #### cross correlate all station pairs ####
    for ii in range(len(st[:-1])):
        for jj in range(ii+1,len(st)):
            index,value,cc_vector=xcorr(st[ii],st[jj],shift_len=cc_shift_length,full_xcorr=True)
            #### if best xcorr value is negative, find the best positive one ####
            if value<0:
                index=cc_vector.argmax()-cc_shift_length
                value=cc_vector.max()
            dt = index/st[0].stats.sampling_rate
            #### check that the best lag is at least the vmin value
            #### and check for minimum cross correlation value
            all_vmin=np.array(vmin).min()
            if np.abs(dt) < intsd[ii,jj]/all_vmin and value > min_cc:
                lags   = np.append(lags,dt)
                lags_inds1 = np.append(lags_inds1,ii)
                lags_inds2 = np.append(lags_inds2,jj)

    #### return lag times, and 
    return lags, lags_inds1, lags_inds2


def associator(lags_inds1,lags_inds2,st, min_chan):
    from itertools import combinations
    #### successively try to associate, starting with all stations
    #### and quit at config.min_sta

    counter   = 0
    mpk = len(st)

    while counter==0 and mpk>=min_chan:
        cmbm = np.array(list(combinations(range(0,len(st)),mpk)))
        cntr = len(cmbm)
        # find how many relevant picks exist for all combinations of delay times
        ncntrm = np.zeros((mpk,mpk,cntr))

        for jj,trig in enumerate(lags_inds1):
            for ii in range(0,cntr):
                if np.sum(lags_inds1[jj] == cmbm[ii,]) == 1 and np.sum(lags_inds2[jj] == cmbm[ii,]) == 1:
                    ind1=(lags_inds1[jj] == cmbm[ii,]).argmax()
                    ind2=(lags_inds2[jj] == cmbm[ii,]).argmax()
                    ncntrm[ind1,ind2,ii] = 1

        # if one of the row/column sums is at least 3, accept it
        cmbm2  = np.zeros((cntr,mpk))
        cmbm2n = np.zeros(cntr)

        for ii in range(cntr):
            if np.sum(np.sum(ncntrm[:,:,ii],1) == 0) == 1:
                cmbm2[counter,:] = cmbm[ii,:]
                # total number of qualifying picks
                cmbm2n[counter] = np.sum(np.sum(ncntrm[:,:,ii],1))
                counter = counter + 1
        # if no matches, decrement and try again
        if counter==0:
            mpk = mpk -1
    
    cmbm2 = cmbm2.astype('int')
    cmbm2n = cmbm2n.astype('int')

    return cmbm2, cmbm2n, counter, mpk


def inversion(cmbm2n,cmbm2,intsd,ints_az,lags_inds1,lags_inds2,lags,mpk):
    # for jj in range(counter):
    jj=0
    # the size of the dt and Dm3
    dt  = np.zeros(cmbm2n[jj])
    Dm3 = np.zeros((cmbm2n[jj],2))

    # initialize interstation distance and azimuth vectors
    ds = np.array([])
    az = np.array([])

    # grab interstation distance and azimuth for all pairs in this tuple
    for num,kk in enumerate(cmbm2[jj,range(0,mpk-1)]):
        for ii in cmbm2[jj,range(num+1,mpk)]:
            ds = np.append(ds,intsd[kk,ii])
            az = np.append(az,ints_az[kk,ii])

    # some counters to find if there is a match in the trgs vector
    mtrxc = 0
    dacnt = 0

    # all 5 may not exist
    for kk in range(0,mpk-1):
        for ii in range(kk+1,mpk):
            tmp=np.array([lags_inds1,lags_inds2]).T - np.repeat(np.array([cmbm2[jj,kk],cmbm2[jj,ii]],ndmin=2),len(lags_inds1),0)
            tmp=np.sum(np.abs(tmp),1)
            mmin = tmp.min()
            mloc = tmp.argmin()
            if mmin==0:
                dt[mtrxc] = lags[mloc]
                Dm3[mtrxc,:] = [ds[dacnt]*np.cos(az[dacnt]*(np.pi/180.0)) , ds[dacnt]*np.sin(az[dacnt]*(np.pi/180.0))]
                mtrxc=mtrxc+1
            dacnt=dacnt+1
    Dm3=Dm3/1000.0  # convert to kilometers

    # generalized inverse of slowness matrix
    Gmi = np.linalg.inv(np.matmul(Dm3.T,Dm3))
    # slowness - least squares
    sv = np.matmul(np.matmul(Gmi,Dm3.T),dt.T)
    # velocity from slowness
    velocity = 1/np.sqrt(np.square(sv[0])+np.square(sv[1]))
    # cosine and sine for backazimuth
    caz3 = velocity*sv[0]
    saz3 = velocity*sv[1]
    # 180 degree resolved backazimuth to source
    azimuth = np.arctan2(saz3,caz3)*(180/np.pi)
    if azimuth<0:
        azimuth=azimuth+360
    # rms
    rms = np.sqrt(np.mean(np.square(np.matmul(Dm3,sv)-dt.T)))

    return velocity, azimuth, rms

    
def xcorr_align_stream(st, cc_shift_length=3*50):

    shift_len=cc_shift_length
    shifts=[]
    for i,tr in enumerate(st):
        a,b,c=xcorr(st[0],tr,shift_len,full_xcorr=True)
        if b<0:
            a=c.argmax()-shift_len
        shifts.append(a/tr.stats.sampling_rate)

    group_streams=Stream()
    T1=st[0].copy().stats.starttime
    T2=st[0].copy().stats.endtime
    for i, tr in enumerate(st):
        tr = tr.copy().trim(tr.stats.starttime-shifts[i],tr.stats.endtime-shifts[i],
            pad=True, fill_value=0)
        tr.trim(tr.stats.starttime+1,tr.stats.endtime-1,pad=True,fill_value=0)
        tr.stats.starttime=T1
        group_streams += tr

    ST=st[0].copy()
    for tr in st[1:]:
        ST.data=ST.data+tr.data
    ST.data=(ST.data/len(st))*config.digouti
    ST.trim(T1,T2)
    return ST