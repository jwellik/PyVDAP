def infra_seis_figure(server, port, st, volcano, T0, config, mx_pressure, infra_duration=600, seis_duration=None):
    import matplotlib as m
    m.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import LinearSegmentedColormap
    from PIL import Image
    import matplotlib.dates as mdates
    
    from pywellik.pyvdap.pywinston import grab_data
    
    seis_duration = infra_duration if seis_duration is None else seis_duration
    start = time.time()
    ##### get seismic data #####
    seis = grab_data(volcano['seismic_scnl'],T0-seis_duration, T0,fill_value='interpolate') # check
    ##### get infrasound data #####
    infra_scnl = ['{}.{}.{}.{}'.format(tr.stats.station,tr.stats.channel,tr.stats.network,tr.stats.location) for tr in st]
    infra = grab_data(infra_scnl,T0-infra_duration, T0,fill_value='interpolate') # check
    end = time.time()
    print('{:.2f} seconds to grab figure data.'.format(end - start))


    ###################################################
    ################# plot infrasound #################

    #### preprocess data ####
    infra.detrend('demean')
    infra.taper(max_percentage=None,max_length=config.taper_val)
    infra.filter('bandpass',freqmin=config.f1,freqmax=config.f2)
    [tr.decimate(2,no_filter=True) for tr in infra if tr.stats.sampling_rate==100]
    [tr.decimate(2,no_filter=True) for tr in infra if tr.stats.sampling_rate==50]
    [tr.resample(25) for tr in infra if tr.stats.sampling_rate!=25]

    ##### stack infrasound data #####
    stack=xcorr_align_stream(infra,config)

    ##### plot stack spectrogram #####
    plt.figure(figsize=(4.5,4.5))
    colors=cm.jet(np.linspace(-1,1.2,256))
    color_map = LinearSegmentedColormap.from_list('Upper Half', colors)

    ax=plt.subplot(len(seis)+3,1,1)
    ax.set_title(config.alarm_name+' Alarm: '+volcano['volcano']+ ' detection!')
    stack.spectrogram(title='',log=False,samp_rate=25,dbscale=True,per_lap=0.7,mult=25.0,wlen=3,cmap=color_map,axes=ax)
    ax.set_yticks([3,6,9,12])
    ax.set_ylim(0,12.5)
    ax.set_ylabel(stack.stats.station+'\nstack',fontsize=5,
                                         rotation='horizontal',
                                         multialignment='center',
                                         horizontalalignment='right',
                                         verticalalignment='center')
    ax.yaxis.set_ticks_position('right')
    ax.tick_params('y',labelsize=4)
    ax.set_xticks([])

    ##### plot stack trace #####
    ax=plt.subplot(len(seis)+3,1,2)
    t1=mdates.date2num(infra[0].stats.starttime.datetime)
    t1=round(t1*24*60)/(24*60) # round to nearest minute
    t2=mdates.date2num(infra[0].stats.endtime.datetime)
    t2=round(t2*24*60)/(24*60) # round to nearest minute
    t_vector=np.linspace(t1,t2,stack.stats.npts)
    plt.plot(t_vector,stack.data,color='k',LineWidth=0.2)
    ax.set_ylabel(stack.stats.station+'\nstack',fontsize=5,
                                         rotation='horizontal',
                                         multialignment='center',
                                         horizontalalignment='right',
                                         verticalalignment='center')
    ax.yaxis.set_ticks_position('right')
    ax.tick_params('y',labelsize=4)
    ax.set_xlim(t1,t2)
    t_ticks=np.linspace(t1,t2,6)
    ax.set_xticks(t_ticks)
    ax.set_xticklabels([mdates.num2date(t).strftime('%H:%M') for t in t_ticks])
    ax.tick_params('x',labelsize=5)
    ax.set_xlabel( '{:.0f} Minute Infrasound Stack\nPeak Pressure: {:.3f} Pa'.format(round((t2-t1)*24*60), mx_pressure) )
    #ax.set_xlabel('{:.0f} Minute Infrasound Stack\n{} UTC,   Peak Pressure: {:.1f} Pa'.format(round((t2-t1)*24*60),
    #                                                    tr.stats.starttime.strftime('%Y-%b-%d'),
    #                                                    mx_pressure))
    ###################################################
    ###################################################


    ###################################################
    ################## plot seismic ###################

    #### preprocess data ####
    seis.detrend('demean')
    [tr.decimate(2,no_filter=True) for tr in seis if tr.stats.sampling_rate==100]
    [tr.decimate(2,no_filter=True) for tr in seis if tr.stats.sampling_rate==50]
    [tr.resample(25) for tr in seis if tr.stats.sampling_rate!=25]


    for i,tr in enumerate(seis):
        ax=plt.subplot(len(seis)+3,1,i+1+3)
        tr.spectrogram(title='',log=False,samp_rate=25,dbscale=True,per_lap=0.5,mult=25.0,wlen=6,cmap=color_map,axes=ax)
        ax.set_yticks([3,6,9,12])
        ax.set_ylabel(tr.stats.station+'\n'+tr.stats.channel,fontsize=5,
                                                             rotation='horizontal',
                                                             multialignment='center',
                                                             horizontalalignment='right',
                                                             verticalalignment='center')
        ax.yaxis.set_ticks_position('right')
        ax.tick_params('y',labelsize=4)

        if i!=len(seis)-1:
            ax.set_xticks([])
        else:
#            d_sec=np.linspace(0,3600,7) # JJW2 commented out
            d_sec=np.linspace(0,seis_duration,6) # JJW2
            ax.set_xticks(d_sec)
            T=[tr.stats.starttime+dt for dt in d_sec]
            ax.set_xticklabels([t.strftime('%H:%M') for t in T])
            ax.tick_params('x',labelsize=5)
            ax.set_xlabel('{:.0f} Minute Local Seismic Data'.format(round(tr.stats.endtime-tr.stats.starttime)/60))

    ###################################################
    ###################################################

    plt.subplots_adjust(left=0.08,right=.94,top=0.92,bottom=0.1,hspace=0.1)
    filename=utils.tmp_figure_dir+'/'+UTCDateTime.utcnow().strftime('%Y%m%d_%H%M%S_%f')
    print(filename)
    plt.savefig(filename+'.png',dpi=250,format='png')
    #m=Image.open(filename) # JJW2
    #remove(filename) # JJW2
    #filename=filename+'.jpg' # JJW2
    #im.save(filename) # JJW2
    print('file saved')

    return filename