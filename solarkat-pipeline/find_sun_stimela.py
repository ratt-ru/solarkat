#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk

import numpy as np
import logging
import numpy
import sys, os
from pyrap.tables import table
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon

from  MSUtils.msutils import addcol
from  MSUtils.msutils import copycol


def deg2str(angle,ra,hour=False):
  if hour:
      angle /= 15
  deg = int(angle)
  angle = (angle - deg)*60*(-1 if ra<0 else 1)
  min = int(angle)
  angle = (angle - min)*60
  return f"{deg}{'h' if hour else 'd'}{min}m{angle:.2f}s"

def backup_UVW(ms,UVW_colname):
    addcol(ms, colname='UVW_backup',clone=UVW_colname)
    copycol(ms,UVW_colname,'UVW_backup')
    print('UVWs saved successfuly.')

def restore_UVW(ms,UVW_colname):
    addcol(ms, colname=UVW_colname,clone='UVW_backup')
    copycol(ms,UVW_colname,'UVW_backup')
    print('UVWs restored successfuly.')

def sun_coordinates(ms,outfile):

    def format_coords(ra0,dec0):
        c = SkyCoord(ra0*u.deg,dec0*u.deg,frame='fk5')
        hms = str(c.ra.to_string(u.hour))
        dms = str(c.dec)
        return hms,dms

    # MeerKAT
    obs_lat = -30.71323598930457
    obs_lon = 21.443001467965008
    loc = EarthLocation.from_geodetic(obs_lat,obs_lon) #,obs_height,ellipsoid)
    maintab = table(ms,ack=False)
    scans = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    lines=[]
    for scan in scans:
        subtab = maintab.query(query='SCAN_NUMBER=='+str(scan))
        t_scan = numpy.mean(subtab.getcol('TIME'))
        t = Time(t_scan/86400.0,format='mjd')
        with solar_system_ephemeris.set('builtin'):
            sun = get_body('Sun', t, loc)
            sun_ra = sun.ra.value
        sun_dec = sun.dec.value
        sun_hms=format_coords(sun_ra,sun_dec)
        print(sun_hms[0],sun_hms[1])
        lines.append(sun_hms[0]+' '+sun_hms[1])
    print(lines)

    with open(outfile,'wt') as f:
        for line in lines:
           f.write(line)
           f.writelines('\n')
 
def shift_to_sun(ms,sun_coords,splitted_ms_dir):
    maintab = table(ms,ack=False)
    scans = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    for scan in scans:
      splitted_ms=splitted_ms_dir+ ms.replace(".ms","_scan"+str(scan)+".ms")
      with open(sun_coords) as f:
          lines = f.readlines()
          for line in lines:
            os.system(f"chgcentre {splitted_ms} {line}")
            print('scan'+str(scan)+ ' Done')
            print(UVW_new,'Old UVW are:')
            UVW_new=maintab.getcol('UVW')
            print(UVW_new,'New UVW are:')
            UVW_old=maintab.getcol('UVW_backup')


def perscan_timerange_intervals(interval,scan_list):

    dic={ "-": "/", " ":"/"}

    def replace_all(text, dic):
        
        for i, j in dic.items():
            text = text.replace(i, j)
        return text

    array=[]
    timerange_array=[]
    tb = table(scan_list)
    all_times = list(np.unique(tb.getcol('TIME')))
    t0 = all_times[0]
    t1 = all_times[-1]
    dt = (t1-t0)/(interval)
    for i in range(interval):
        t2=dt*i+t0
        t_iso = Time(t2/86400.0,format='mjd').iso
        array.append(t_iso)
    for i in range(len(array)):
        if i < (len(array)-1):
            timerange=replace_all(array[i],dic)+'~'+replace_all(array[i+1],dic)
            timerange_array.append(timerange)
        else:
            print()
    timeranges=timerange_array
    print(timerange_array)
            
  
            
   
