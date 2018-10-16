## netstat.netstat.__init__

def intro_message():
    import os
    print("NETSTAT ({})".format(os.getcwd()))
    
    
def print_network_info(config):
    print("Network name : {}".format(config['DEFAULT']['name']))
    print(" Minimum gap duration      : {} sec".format(config['thresholds']['min_gap_dur']))
    print(" Tranmission Acceptability : {}%".format(config['thresholds']['sta_acceptability']))
    print(" Minimum sations           : {}".format(config['thresholds']['min_sta']))
    print("")
    
def print_station_status(line, sta_status_str, ngaps, total_gap_duration, data_percent):
    print(line)
    print('STATUS = {}'.format(sta_status_str))
    print('Gaps      : {}'.format(ngaps))
    print('Gap dur   : {:d} sec'.format(int(total_gap_duration)))
    print('Data %    : {:4.1f}%'.format(data_percent))
    print(' ')
    
def print_network_status(name, network_status_str, good_sta, min_sta, total_sta, network_status_msg):
    print('-'*45)
    print('{} network status: {} ({} of {})'.format(name, network_status_str, good_sta, total_sta))
    print('   (BAD < {} <= ADEQUATE < {} <= GOOD)'.format( int(min_sta[0]), int(min_sta[1]) ))
    print('-'*45)
    [print(m) for m in network_status_msg]