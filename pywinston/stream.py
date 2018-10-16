### Based on miscellaneous codes by Aaron Wech
from obspy import Stream
import numpy as np


def remove_empty_tr(st, min_tr=None):
    '''Removes empty traces
    If min_tr is specified, returns a boolean of whether or not there are enough populated
     traces in the stream
    
    # Remove empty traces and return the stream object
    >>> st = remove_empty_tr( st )
    
    # Removes empty traces and returns boolean if there are at least 3 traces remaining
    >>> remove_empty_tr( st, min_tr = 3 )
    False
    '''
    
    for tr in st:
        if np.sum(np.abs(np.abs(tr.data)))==0:
            st.remove(tr)

    if min_tr is not None:
            st = True if len(st)>=min_tr else False
            
    return st


def remove_gappy_tr(st):
    '''Removes gappy traces
    '''

    for tr in st:
        num_zeros=len(np.where(tr.data==0)[0])
        if num_zeros/float(tr.stats.npts)>0.01:
            st.remove(tr)
            
    return st


def preprocess(st, taper_val=5.0, bpass=(0.4, 10)):
    '''Routine preprocess of stream data
    Demeans, tapers, bandpass filters, and downsamples data
    
    Parameters
    ----------
    st          : ObsPy Stream object
    taper_val   : int
        Passed to Stream.taper as the 'max_length' argument
    bpass       : tuple
        low pass and hi pass corner for the bandpass filter (e.g., (1.0, 10.0) )
    '''
    #### preprocess data ####
    st.detrend('demean')
    st.taper(max_percentage=None,max_length=taper_val)
    st.filter('bandpass',freqmin=bpass[1],freqmax=bpass[1])
    for tr in st:
        if tr.stats['sampling_rate']==100:
            tr.decimate(2)
        if tr.stats['sampling_rate']!=50:
            tr.resample(50.0)
    
    return st


def check_amplitude_thresh( st, digouti, min_pa, min_cha ):
    st=Stream([tr for tr in st if np.any(np.abs(tr.data*digouti)>min_pa)])
    if len(st)<min_chan:
        return False
    else:
        return True