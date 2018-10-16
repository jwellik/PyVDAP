def craft_and_send_email(t1,t2,config,volcano,d_Azimuth,velocity,mx_pressure,filename):
    from pandas import Timestamp
    # create the subject line
    subject='{} Airwave Detection'.format(volcano['volcano'])

    # create the test for the message you want to send
    message='{} alarm:\n'.format(config.alarm_name)
    message='{}{} detection!\n\n'.format(message,volcano['volcano'])
    message='{}Start: {} (UTC)\nEnd: {} (UTC)\n\n'.format(message,t1.strftime('%Y-%m-%d %H:%M'),t2.strftime('%Y-%m-%d %H:%M'))
    t1_local=Timestamp(t1.datetime,tz='UTC')
    t2_local=Timestamp(t2.datetime,tz='UTC')
    t1_local=t1_local.tz_convert('US/Alaska')
    t2_local=t2_local.tz_convert('US/Alaska')
    message='{}Start: {} ({})'.format(message,t1_local.strftime('%Y-%m-%d %H:%M'),t1_local.tzname())
    message='{}\nEnd: {} ({})\n\n'.format(message,t2_local.strftime('%Y-%m-%d %H:%M'),t2_local.tzname())

    message='{}d_Azimuth: {:+.1f} degrees\n'.format(message,d_Azimuth)
    message='{}Velocity: {:.0f} m/s\n'.format(message,velocity*1000)
    message='{}Max Pressure: {:.1f} Pa'.format(message,mx_pressure)

    utils.send_alert(config.alarm_name,subject,message,filename)
    #utils.post_mattermost(subject,message,config.alarm_name,filename) #JJW2
    # delete the file you just sent
    if filename:
        remove(filename)
