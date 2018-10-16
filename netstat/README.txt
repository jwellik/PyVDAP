IDEAL DESCRIPTION OF FUNCTIONALITY
(This is what the documentation would look like if the program did everything exactly
the way I wante it to).

DIRECTORY STRUCTURE

./netstat
    /configs
        Agung.config
    /reports
        Agung.txt
        Agung.html
        Agung.png
    /tmp                   # temporary location of txt file w data
    netstat.py
    netstat.sh
    plot_netstatus.py

USAGE
"Alert mode" Run NETALERT for a configuration file.
This is the best mode for regular updates to send e-mail alerts and check the recent status
$ python netalert.py Agung

"History mode" Run NETALERT for a configuration file with a start and stop time.
This over-rides [DEFAULT]['duration'] and instead produces a report over the given time range.
E-mail alerts and status updates to the configuration file are *not* produced.
Start/Stop dates must be 'YYYY-mm-dd'. 'now' may be used as the stop date.
$ python netalert.py Agung 2017-10-01 now
