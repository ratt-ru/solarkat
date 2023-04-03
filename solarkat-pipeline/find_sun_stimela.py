#!/usr/bin/env python
import re
import glob
import numpy
import logging
import sys, os
import subprocess
import numpy as np
from astropy.time import Time
from pyrap.tables import table
from astropy import units as u
from casacore.tables import table
from  MSUtils.msutils import addcol
from  MSUtils.msutils import copycol
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon

 

def hms2deg(hms):
    """
    Function to convert hms to degrees
    """
    hms_angle = Angle(hms, unit='hour')
    return hms_angle.degree


def dms2deg(dms):
    """
    Function to convert dms to degrees
    """
    dms_angle = Angle(dms, unit='degree')
    return dms_angle.degree


def scan_numbers(ms):
    """
    Extracts unique scan numbers from the Measurement Set
    """
    from casacore.tables import table
    import numpy
    scans=[]
    maintab = table(ms,ack=False)
    scan_no = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    for scan in scan_no:
        scans.append(str(scan))
    maintab.close()
    return scans

#______________________________________________________________________________________________________________________________________________

def perscan_timerange_intervals1(interval,scan_list):
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
    tb.close()
    return timeranges
    
    return timeranges


def perscan_timerange_intervals(interval, scan_list):
    array = []
    timerange_array = []
    tb = table(scan_list)
    all_times = list(np.unique(tb.getcol('TIME')))
    t0 = all_times[0]
    t1 = all_times[-1]
    array = [Time(t/86400.0,format='mjd').iso.replace(" ", "/").replace("-", "/") for t in np.linspace(t0, t1, interval + 1)]
    for i in range(len(array) - 1):
        timerange = array[i] + "~" + array[i+1]
        timerange_array.append(timerange)
    timeranges = timerange_array
    tb.close()
    return timeranges


#______________________________________________________________________________________________________________________________________________

def get_old_coords(ms_list, outfile):
    """
    Extracts the old coordinates (ra and dec) of the foV from a measurement set using the PHASE_DIR column and writes them to a file.

    Parameters:
    ms_list (str): Path to the measurement sets files.
    outfile (str): Path to the output file.

    Returns:
    None
    """

    old_coords = []
    for ms in ms_list:

        # 1. Read the PHASE_DIR column, the values in radians
        field_dir = table(f"{ms}::FIELD", readonly=True)
        phase_dir = field_dir.getcol("PHASE_DIR")

        # Extract the RA and Dec values from the PHASE_DIR column
        ra = phase_dir[0,0,0]
        dec = phase_dir[0,0,1]

        #2.convert to sky coordinates hms/dms
        sk = SkyCoord(ra*u.rad, dec*u.rad)

        # format th ra to hms
        ra_hms = sk.ra.to_string(unit=u.hourangle, pad=True)#, precision=1)

        # format dec to dms
        dec_dms = sk.dec.to_string(unit=u.degree, pad=True, alwayssign= True)#,precision=1)
        old_coords.append(f"{ra_hms} {dec_dms}")
        print(old_coords)
        field_dir.close()

    with open(outfile,'wt') as f:
        for old_coord in old_coords:
            f.write(old_coord)
            f.writelines('\n')

    return old_coords

#_______________________________________________________________________________________________________________________________________________


def get_sun_coordinates(ms, outfile):
    """
    Extracts the coordinates of the Sun from a measurement set and writes them to a file.

    Parameters:
    ms (str): Path to the measurement set file.
    outfile (str): Path to the output file.

    Returns:
    None
    """
    def format_coords(ra0,dec0):
        c = SkyCoord(ra0*u.deg,dec0*u.deg,frame='fk5')
        hms = str(c.ra.to_string(u.hour))
        dms = str(c.dec)
        return hms,dms

    # MeerKAT
    obs_lat = -30.71323598930457
    obs_lon = 21.443001467965008
    loc = EarthLocation.from_geodetic(obs_lat, obs_lon)

    maintab = table(ms, ack=False)
    scans = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    lines = []

    for scan in scans:
        subtab = maintab.query(query='SCAN_NUMBER==' + str(scan))
        t_scan = numpy.mean(subtab.getcol('TIME'))
        t = Time(t_scan / 86400.0, format='mjd')

        with solar_system_ephemeris.set('builtin'):
            sun = get_body('Sun', t, loc)
            sun_ra = sun.ra.value
            sun_dec = sun.dec.value
            sun_hms, sun_dms = format_coords(sun_ra, sun_dec)
            lines.append(f"{sun_hms} {sun_dms}")

    maintab.close()

    with open(outfile, 'wt') as f:
        for line in lines:
            f.write(line + '\n')

#____________________________________________________________________________________________________________________________________________


def get_pertime_sun_coordinates(ms_list, outfile):

    """
    Extracts the scan coordinates of the Sun from a list of  measurement sets splitted by time intervals and writes them to a file.

    Parameters:
    ms_list (str): Path to the measurement sets pertime files.
    outfile (str): Path to the output file.

    Returns:
    None
    """

    def format_coords(ra0,dec0):
        c = SkyCoord(ra0*u.deg,dec0*u.deg,frame='fk5')
        hms = str(c.ra.to_string(u.hour))
        dms = str(c.dec)
        return hms,dms

    # MeerKAT
    obs_lat = -30.71323598930457
    obs_lon = 21.443001467965008
    loc = EarthLocation.from_geodetic(obs_lat,obs_lon)
    lines=[]
    for ms in ms_list:
        maintab = table(ms, ack=False)
        scans = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
        for scan in scans:
            subtab = maintab.query(query='SCAN_NUMBER=='+str(scan))
            t_scan = numpy.mean(subtab.getcol('TIME'))
            t = Time(t_scan/86400.0,format='mjd')
            with solar_system_ephemeris.set('builtin'):
                sun = get_body('Sun', t, loc)
                sun_ra = sun.ra.value
                sun_dec = sun.dec.value
                sun_hms,sun_dms=format_coords(sun_ra,sun_dec)
                lines.append(f"{sun_hms} {sun_dms}")
        maintab.close()

        with open(outfile, 'wt') as f:
            for line in lines:
                f.write(line + '\n')
                
#___________________________________________________________________________________________________________________________________________


def rename_model_data_column(ms_list):
    """
    Rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL' in a list of MS tables.
    """
    for ms_name in ms_list:
        # open the MS table in read-write mode
        ms = table(ms_name, readonly=False)
        # rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL'
        ms.renamecol('MODEL_DATA', 'MODEL_DATA_ORIGINAL')
        # close the MS table
        ms.close()

#____________________________________________________________________________________________________________________________________________


def shift_coordinates(ms_list, coords, splitted_ms_dir, datacolumn):
    """
    A funtion that takes a list of scans and coordinates shift/rephase it for a specific colunm (CORRECTED_DATA) and iutput in the splitted_ms_dir directory 
    Parameters:
    ms_list (str): Path to the measurement sets files.
    coords (File): Path to the coordinate files.
    splitted_ms_dir (Directory): Path to the directory
    datacolumn (str): Datacolumn to use
    """
    with open(coords) as f:
        lines = f.readlines()
        num_lines = len(lines)
    for i, ms in enumerate(ms_list):
        splitted_ms = os.path.join(splitted_ms_dir, os.path.basename(ms).replace(".ms", ".ms"))
        with open(coords) as f:
            line = lines[i % num_lines].strip()
            os.system(f"chgcentre {ms} {line} {splitted_ms} datacolumn={datacolumn}")
        print(f'Scan {ms} Done.')

#___________________________________________________________________________________________________________________________________________


def add_column_to_ms(ms_list, col_name, like_col, like_type=None):
    """
    Add a column to each measurement file in the list of measurement files
    """
    success = False
    for ms in ms_list:
        try:
            tb = casacore.tables.table(ms)
        except Exception:
            continue
        if col_name not in tb.colnames():
            """
            Get column description from column 'like_col'
            """
            desc = tb.getcoldesc(like_col)

            desc[str('name')] = str(col_name)
            desc[str('comment')] = str(desc['comment'].replace(" ", "_"))  # got this from Cyril, not sure why
            dminfo = tb.getdminfo(like_col)
            dminfo[str("NAME")] =  "{}-{}".format(dminfo["NAME"], col_name)
            # if a different type is specified, insert that
            if like_type:
                desc[str('valueType')] = like_type
            tb.addcols(desc, dminfo)
            tb.close()
            success = True
    return success

#____________________________________________________________________________________________________________________________________________



def rename_columns(ms_list):
    """
    Rename columns in a list of MS tables.
    """
    for ms_name in ms_list:
        """
        open the MS table in read-write mode
        """
        ms = table(ms_name, readonly=False)
        """
        rename the 'MODEL_DATA' column to 'SUN_MODEL'
        """
        if 'MODEL_DATA' in ms.colnames():
            ms.renamecol('MODEL_DATA', 'SUN_MODEL')
        # rename the 'MODEL_DATA_ORIGINAL' column to 'MODEL_DATA'
        if 'MODEL_DATA_ORIGINAL' in ms.colnames():
            ms.renamecol('MODEL_DATA_ORIGINAL', 'MODEL_DATA')
        # close the MS table
        ms.close()

#____________________________________________________________________________________________________________________________________________


def create_ds9_region_from_file(input_file, output_dir):
    """
    Function to create a ds9 region file for each coordinate in the input file
    """
    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            ra, dec = map(str.strip, line.split())
            ra = hms2deg(ra)
            dec = dms2deg(dec)
            create_ds9_region(ra, dec, i, output_dir)

def create_ds9_region(ra, dec, index, output_dir):
    # Function to create a ds9 region file
    sun_region = f'''# Region file format: DS9 CARTA 3.0.0-beta.3 
global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle({ra}, {dec}, 1188.0000") # color=green
'''

    with open(f'{output_dir}/1671435077_sdp_l0_1024ch_GRS1747-312_scan_{index}-MFS-image.reg', 'w') as f:
        f.write(sun_region)






