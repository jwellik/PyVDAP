import os

import configparser

import requests
import urllib.request
from urllib.parse import urlparse, parse_qs, parse_qsl, urlencode

import pandas as pd
import numpy as np

import warnings

import netstat

##############################
###
##############################

    ### read config file ###
config_name = 'Agung' # first input argument

config_pathname = '/configs/{}.config'.format(config_name)
config = configparser.ConfigParser(delimiters=('='))
config.read(os.getcwd() + config_pathname)


    ### read network and server information ###
name = config['DEFAULT']['name']
server = config['DEFAULT']['server']
port = config['DEFAULT']['port']
duration = config['DEFAULT']['duration']
outdir = config['DEFAULT']['outdir']
print('{} @ {}:{}'.format(name,server,port))


    ### read threshold information ###
sta_acceptability = config['thresholds']['sta_acceptability']
min_sta = config['thresholds']['min_sta']
mgd = config['thresholds']['min_gap_dur']

prev_network_status = config['status']['current_status']
print(sta_acceptability)
print(min_sta)

stations = config['stations']['scnl'].splitlines()
stations

network_status = []
network_status_msg = []
for line in stations:
    if line is not '':
        line = line.replace(' ','').split(',')
        scnl = line[0]
        
            # initialize local version of thresholds for each indiviudal station
        sta_acceptability = config['thresholds']['sta_acceptability']
        min_sta = config['thresholds']['min_sta']
        mgd = config['thresholds']['min_gap_dur']
       
            # process station-specific parameters
        for e in line:
            if 'sta_acceptability' in e:
                elem = e.split(':')
                sta_acceptability = elem[1]
        
            ### get gap data from Winston ###
            # geturl example: 'http://vdap.org:16024/gaps?code=TMKS_EHZ_VG_00&t1=-24&wc=1'
        geturl = 'http://{}:{}/gaps?code={}&t1={}&mgd={}&wc=1'.format(server,port,scnl,duration,mgd)
        tmpdata = os.getcwd() + '/tmp/{}.txt'.format(scnl)
        print('Requesting data: {}'.format(geturl))
        try:
            urllib.request.urlretrieve(geturl, tmpdata)
        except:
            warnings.warn('Failed to retrieve data.')
            tmpdata = os.getcwd() + '/tmp/TMKS_EHZ_VG_00_test.txt'
            ### / get gap data from Winston ###
            
            # process station gap data            
        df = pd.read_csv(tmpdata,
                         comment='#', header=None, delim_whitespace=True,
                         names=['GapStart', 'GapEnd', 'Duration'])
        #os.remove(tmpdata)
        ngaps = df.shape[0]
        total_gap_duration = df['Duration'].sum()
        total_dur = 86400
        gap_percent = total_gap_duration / total_dur * 100
        data_percent = 100 - gap_percent
        sta_status_str = 'GOOD' if data_percent >= float(sta_acceptability) else 'BAD'
        sta_status = 1 if sta_status_str is 'GOOD' else 0
        network_status.append(sta_status)
        network_status_msg.append('{} : {:4.1f}% ({:>4})'.format(
            scnl.upper() if sta_status else scnl.lower(),
            data_percent,
            sta_status_str))
        
            # print station message
        print(line)
        print('STATUS = {}'.format(sta_status_str))
        print('Gaps      : {}'.format(ngaps))
        print('Gap dur   : {:d} sec'.format(int(total_gap_duration)))
        print('Data %    : {:4.1f}%'.format(data_percent))
        print(' ')

    # process network results
good_sta = int(np.array(network_status).sum())
total_sta = len(stations) - 1 # makes up for stations[0] being blank
network_status_str = 'GOOD' if good_sta >= int(min_sta) else 'BAD'

    # print network results
print('-'*45)
print('{} network status: {} ({} out of {}/{})'.format(name, network_status_str, good_sta, min_sta, total_sta))
print('-'*45)
[print(m) for m in network_status_msg]