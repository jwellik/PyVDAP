import os

import configparser
import argparse
import requests
import urllib.request
from urllib.parse import urlparse, parse_qs, parse_qsl, urlencode
import warnings

import pandas as pd
import numpy as np

import netstat
from netstat.plotting import netstatmap

##############################
### READ INPUT ARGUMENTS
##############################

netstat.intro_message()


parser = argparse.ArgumentParser()
parser.add_argument("config_name", type=str,
                    help="Configuration file name (no extension)")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="produce verbose output")
args = parser.parse_args()

config_pathname = '{}/configs/{}.config'.format(os.getcwd(), args.config_name)
config = configparser.ConfigParser(delimiters=('='))
config.read(config_pathname)

print("Configuration file: {}".format(config_pathname))

##############################
### read network and server information ###
##############################

name = config['DEFAULT']['name']
server = config['DEFAULT']['server']
port = config['DEFAULT']['port']
duration = config['DEFAULT']['duration']
outdir = config['DEFAULT']['outdir']
print('{} @ {}:{}'.format(name,server,port))

    ### read threshold information ###
sta_acceptability = config['thresholds']['sta_acceptability']
min_sta = config['thresholds']['min_sta']
min_sta = min_sta.replace(' ', '').split(',')
mgd = config['thresholds']['min_gap_dur']

prev_network_status = config['status']['current_status']
if args.verbose: netstat.print_network_info(config)

stations = config['stations']['scnl'].splitlines()
stations

network_status = []
network_status_msg = []
for line in stations:
    if line is not '':
        line = line.replace(' ','').split(',')
        scnl = line[0]
        
            # initialize default thresholds for station transmission status
        sta_acceptability = config['thresholds']['sta_acceptability']
        min_sta = config['thresholds']['min_sta']
        min_sta = min_sta.replace(' ', '').split(',')
        mgd = config['thresholds']['min_gap_dur']
       
            # process station-specific parameters
        for e in line:
            if 'sta_acceptability' in e:
                elem = e.split(':')
                sta_acceptability = elem[1]
        
        tmpdata =  os.getcwd() + '/tmp/{}.txt'.format(scnl)
        # tmpdata, gap_data = pywellik.pyvdap.pywinston.gaps.get_gaps(server, port, scnl, t1=duration, mgd=mgd, wc=1, outfile=tmpdata)
            #################################
            ### get gap data from Winston ###
            # geturl example: 'http://vdap.org:16024/gaps?code=TMKS_EHZ_VG_00&t1=-24&wc=1'
        geturl = 'http://{}:{}/gaps?code={}&t1={}&mgd={}&wc=1'.format(server,port,scnl,duration,mgd)
        print('Requesting data: {}'.format(geturl))
        try:
            urllib.request.urlretrieve(geturl, tmpdata)
        except:
            print('Failed to retrieve data. Using stub data')
            warnings.warn('Failed to retrieve data.')
            tmpdata = os.getcwd() + '/tmp/TMKS_EHZ_VG_00_test.txt'
            ### / get gap data from Winston ###
            
            # process station gap data            
        df = pd.read_csv(tmpdata,
                         comment='#', header=None, delim_whitespace=True,
                         names=['GapStart', 'GapEnd', 'Duration'])
        
        ngaps = df.shape[0]
        total_gap_duration = df['Duration'].sum()
        total_dur = 86400
        gap_percent = total_gap_duration / total_dur * 100
        data_percent = 100 - gap_percent
            ### get gap data from Winston ###
            ####/////////////////////////####
        #os.remove(tmpdata)
        
            ### these need to change if pywinston.gaps.get_gaps is enabled (e.g., "data_percent" --> "gap_data['data_percent']")
        sta_status_str = 'GOOD' if data_percent >= float(sta_acceptability) else 'BAD'
        sta_status = 1 if sta_status_str is 'GOOD' else 0
        network_status.append(sta_status)
        network_status_msg.append('{} : {:4.1f}% ({:>4})'.format(
            scnl.upper() if sta_status else scnl.lower(),
            data_percent,
            sta_status_str))
        
            # print station message
        if args.verbose: netstat.print_station_status(line, sta_status_str, ngaps, total_gap_duration, data_percent)


    # process network results
good_sta = int(np.array(network_status).sum())
total_sta = len(stations) - 1 # makes up for stations[0] being blank
if int(min_sta[0]) > good_sta:
    network_status_str = 'BAD'
elif int(min_sta[0]) <= good_sta < int(min_sta[1]):
    network_status_str = 'ADEQUATE'
else:
    network_status_str = 'GOOD'
      

    # print network results
if args.verbose: netstat.print_network_status(name, network_status_str, good_sta, min_sta, total_sta, network_status_msg)
netstatmap(config)