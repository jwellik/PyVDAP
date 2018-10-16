def netstatmap(config, map_type='terrain', zoom=11, scale_control=True):
    '''Produces map of network status
    *** MANY THINGS STILL HARD-CODED FOR FUEGO ***
    '''
    from bokeh.io import output_file, show
    from bokeh.models import ColumnDataSource, GMapOptions, HoverTool
    from bokeh.plotting import gmap

    print("Producing NETSTAT Test Map")

    output_file("./reports/test_netstat_gmap_{}.html".format(config['DEFAULT']['name']))

    map_options = GMapOptions(lat=float(config['DEFAULT']['vlat']), lng=float(config['DEFAULT']['vlon']),
        map_type=map_type, zoom=zoom, scale_control=scale_control)

    tools = 'pan,wheel_zoom,box_select,box_zoom,lasso_select,reset'.split(',')
    hover = HoverTool(tooltips=[
        ("Name", "@scnl"),
        ("Latitude", "@lat"),
        ("Longitude", "@lon"),
    ])
    tools.append(hover)

    # For GMaps to function, Google requires you obtain and enable an API key:
    #
    #     https://developers.google.com/maps/documentation/javascript/get-api-key
    #
    # Replace the value below with your personal API key:
    p = gmap("AIzaSyA4H87CU8D6-ixMwOt4yL0wFk9ojs19Jvw", map_options, title=config['DEFAULT']['name'],
        tools=tools)

    source = ColumnDataSource(
        data=dict(scnl=['FG3', 'FG8', 'FG9'],
              lat=[ 14.44783,  14.5056,  14.4325],
              lon=[-90.842, -90.9468, -90.9359],
              clr=['green', 'red', 'yellow'])
        )

    p.circle(x="lon", y="lat", size=15, fill_color="clr", fill_alpha=0.8, source=source, legend='seismic')
    p.legend.click_policy="hide" # "mute" --> alternative choice

    show(p)