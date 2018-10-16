def read_latlonconfig( file ):
    '''
    Reads a Winston/Swarm readable latlon.config file
    
    Parameters
    ----------
    file        : str
        Winston/Swarm readable latlon.config file

    Returns
    ----------
    stations    : DataFrame
        DataFrame with SCNL, TimeZone, Latitude, Longitue, Height columns
        * SCNL is '.' delimitted
    '''
    
    '''
    
    Developers Note:
    
    Example line from file:
    BTR EHZ VG 00 = Longitude: 115.37636; Latitude:-8.24523; Height: 1436; TimeZone: Asia/Singapore
                  ^          ^          ^
                  |          |          |-- splits each keyword, value pair from other metadata
                  |          |-- splits key from value
                  |-- splits scnl from scnl metadata
                  
                  
    '''

    import os
    import pandas as pd
    
    stations = []
    with open(file) as myfile:
        for line in myfile:
            if len(line.strip()) > 1: # only parse lines that have data
                line.replace('\n', '')
                scnl, metadata = line.partition("=")[::2] # separate scnl from metadata
                d = {}
                d['SCNL'] = scnl.strip()
                d['SCNL'] = d['SCNL'].replace(' ', '.')
                mdlist = metadata.split(';') # separate up each section of metadata
                for md in mdlist:
                    args = md.split(':') # separate key from value
                    d[args[0].strip()] = args[1].strip()
                stations.append(d)
        
    stations = pd.DataFrame(stations)
    return stations



def grab_data2(client, scnl,T1,T2,fill_value=0):
    '''More robust wrapper for ObsPy.Client.get_waveforms()
    
    Parameters
    ----------
    client          : ObsPy Client
        client server & port from which to grab data
    
    scnl
        list of station names (eg. ['PS4A.EHZ.AV.--','PVV.EHZ.AV.--','PS1A.EHZ.AV.--'])
    
    T1, T2          : ObsPy UTCDateTimes
        start, stop of the data
    
    (fill_value)
        (0) | 'latest' | 'interpolate'
        
    Returns
    ----------
    st              : ObsPy Stream object
        Stream object with gaps accounted for
    '''
    
    from obspy import Stream, Trace
    from obspy.clients.earthworm import Client

    print('{} - {}'.format(T1.strftime('%Y.%m.%d %H:%M:%S'),T2.strftime('%Y.%m.%d %H:%M:%S')))
    print('Grabbing data...')

    st=Stream()

    for sta in scnl:
        try:
            tr=client.get_waveforms(sta.split('.')[2], sta.split('.')[0],sta.split('.')[3],sta.split('.')[1], T1, T2,       cleanup=True)
            if len(tr)>1:
                if fill_value==0 or fill_value==None:
                    tr.detrend('demean')
                    tr.taper(max_percentage=0.01)
                for sub_trace in tr:
                    # deal with error when sub-traces have different dtypes
                    if sub_trace.data.dtype.name != 'int32':
                        sub_trace.data=sub_trace.data.astype('int32')
                    if sub_trace.data.dtype!=dtype('int32'):
                        sub_trace.data=sub_trace.data.astype('int32')
                    # deal with rare error when sub-traces have different sample rates
                    if sub_trace.stats.sampling_rate!=round(sub_trace.stats.sampling_rate):
                        sub_trace.stats.sampling_rate=round(sub_trace.stats.sampling_rate)
                print('Merging gappy data...')
                tr.merge(fill_value=fill_value)
        except:
            tr=Stream()
        # if no data, create a blank trace for that channel
        if not tr:
            from obspy import Trace
            from numpy import zeros
            tr=Trace()
            tr.stats['station']=sta.split('.')[0]
            tr.stats['channel']=sta.split('.')[1]
            tr.stats['network']=sta.split('.')[2]
            tr.stats['location']=sta.split('.')[3]
            tr.stats['sampling_rate']=100
            tr.stats['starttime']=T1
            tr.data=zeros(int((T2-T1)*tr.stats['sampling_rate']),dtype='int32')
        st+=tr
    print('Detrending data...')
    st.detrend('demean')
    st.trim(T1,T2,pad=0)
    return st


def grab_data(server, port, scnl,T1,T2,fill_value=0):
    from obspy import Stream
    from obspy.clients.earthworm import Client

    st=Stream()
    client = Client(server, port)
    for sta in scnl:
        station = sta.split('.')[0]
        channel = sta.split('.')[1]
        network = sta.split('.')[2]
        location = sta.split('.')[3]
        print(station, channel, network, location, T1, T2)
        try:
            tr=client.get_waveforms(network, station, location, channel, T1, T2)
            if len(tr)==0:
                tr=create_trace(sta, T1, T2)
            else:
                if len(tr)>1:
                    if fill_value==0 or fill_value==None:
                        tr.detrend('demean')
                        tr.taper(max_percentage=0.01)
                    tr.merge(fill_value=fill_value)
                #tr.trim(T1,T2,pad=0)
                #tr.detrend('demean')
        except Exception as err:
            print(err)
            print("No data found for "+sta)
            tr=create_trace(sta, T1, T2)
        st+=tr
    print(st)
    return st

def create_trace(sta, T1, T2):
    from obspy import Trace
    from numpy import zeros
    tr=Trace()
    tr.stats['station']=sta.split('.')[0]
    tr.stats['channel']=sta.split('.')[1]
    tr.stats['network']=sta.split('.')[2]
    tr.stats['location']=sta.split('.')[3]
    tr.stats['sampling_rate']=100
    tr.stats['starttime']=T1
    tr.data=zeros(int((T2-T1)*tr.stats['sampling_rate']))
    return tr




def export_sac(server, port, scnl, t1, t2, td=24, fill_value=0, scnl_alias=None, td_offset=None, outdir=None):
    '''Exports Winston data as a SAC file
    
    Parameters
    ----------
    client : ObsPy Client
        datasource (eg. Client('volcano.observatory.org', 16022) )
    
    scnl
        list of station names (eg. ['PS4A.EHZ.AV.--','PVV.EHZ.AV.--','PS1A.EHZ.AV.--'])
    
    t1, t2 : ObsPy UTCDateTime
    
    td : int
        If td is specified, splits data up into 'td' hour chunks
        (default=24)
    
    (fill_value)
        (0) | 'latest' | 'interpolate'
        
    (scnl_alias)
        List of aliases for the stations (eg. ['AAAA.EHZ.00.--','BBBB.EHZ.00.--','CCCC.EHZ.00.--'])
        Useful for renaming data for training purposes
        * Must be same size as 'scnl'
    
    (td_offset) : int
        * THIS FEATURE IS NOT YET AVAILABLE *
        Number of days to offset the exported data date
        Useful for editing data for training purposes
    
    (outdir) : str
        directory where files are stored, creates path if it does not exist
        (pwd)
    
    Yields
    ---------
    file
        format -> S_C_N_L_YYYYmmddHHMMSS.sac
    
    '''
    
    import os

    import pandas as pd
    import datetime as dt
    import numpy as np

    import obspy
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream, Trace

    from pywellik.pyvdap.pywinston import grab_data
    
    #####################
    ### Manage Paths
    #####################
    
    if outdir is None:
        outdir = os.getcwd()
    elif not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    ###################
    ### grab data and export
    ###################
    
    for i, sta in enumerate(scnl):
        
        print(sta)
        
        t = t1
        while t < t2:
            
            #print('Retrieving data : {} - {}'.format(sta, t))
        
            tmp = []
            tmp = grab_data(server, port, [sta], UTCDateTime(t), UTCDateTime(t)+td*60*60, fill_value=fill_value)
            tmp = tmp[0]
            print(tmp)
            if scnl_alias is not None:
                tmp.stats.station = scnl_alias[i].split('.')[0]
                tmp.stats.channel = scnl_alias[i].split('.')[1]
                tmp.stats.network = scnl_alias[i].split('.')[2]
                tmp.stats.location = scnl_alias[i].split('.')[3]
                
            if td_offset is not None:
                print('This feature is not available yet.')
    
            file = '{}_{}_{}_{}-{}'.format(tmp.stats.station,tmp.stats.channel,tmp.stats.network,tmp.stats.location, tmp.stats.starttime.strftime('%Y%m%d%H%M%S'))
            file = os.path.join(outdir,file)
            tmp.write(file, format='sac')
            print('Retrieved : {} {} | Exported : {}'.format(sta, t, file))
            
            t = t + td * 60 * 60
