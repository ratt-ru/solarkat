#!/usr/bin/env python
import re
import numpy
import sys, os
import subprocess
import numpy as np
from astropy.time import Time
from pyrap.tables import table
from astropy import units as u
from  MSUtils.msutils import addcol
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon


#________________________________________________________________________________________________________________________________________________________

#Auxiliar functions

def hms2deg(hms):
    '''
    Function to convert hms to degrees
    '''
    hms_angle = Angle(hms, unit='hour')
    return hms_angle.degree


def dms2deg(dms):
    '''
    Function to convert dms to degrees
    '''
    dms_angle = Angle(dms, unit='degree')
    return dms_angle.degree


def extract_scan_number(ms_scan):
    # Extract the scan number from the scan file name using re
    scan_number = re.search(r"scan_(\d+)\.ms", ms_scan).group(1)
    return int(scan_number)

#________________________________________________________________________________________________________________________________________________________


def extract_and_save_scan_numbers(ms, output_file):
    '''
    Extracts unique scan numbers from the Measurement Set
    measurement_set (str): Measurement Set path
    output_file (str): Output file path to save the scan numbers
    Returns:
    List[str]: List of unique scan numbers as strings
    '''
    print("Extracting scan numbers from {}...".format(ms))

    #print(f"Extracting scan numbers from {ms} ...")
    
    scans = []
    
    with table(ms, ack=False) as maintab:
        scan_numbers = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
        scans = [str(scan) for scan in scan_numbers]

    print("Unique scan numbers extracted:")
    print(scans)

    # Save scan numbers to the output file
    with open(output_file, 'w') as file:
        file.write('\n'.join(scans))

    print("Scan numbers saved to {}...".format(output_file))
    
    return scans

#________________________________________________________________________________________________________________________________________________________

def get_old_coords(ms_list, output_file):
    '''
    Extracts the old coordinates (ra and dec) of the foV from a measurement set using the PHASE_DIR column and writes them to a file.

    Parameters:
    ms_list (str): Path to the measurement sets files.
    outfile (str): Path to the output file.

    Returns:
    None
    '''
    sorted_ms_list = sorted(ms_list, key=lambda x: int(x.split('_scan_')[1].split('.')[0]))

    old_coords = []
    for ms in sorted_ms_list:

        print("Processing {}...".format(ms))
        # 1. Read the PHASE_DIR column, the values in radians
        field_dir = table("{}::FIELD".format(ms), readonly=True)
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
        old_coords.append("{} {}".format(ra_hms, dec_dms))
        field_dir.close()
        print("Processing of {} complete.".format(ms))

    with open(output_file,'wt') as f:
        for old_coord in old_coords:
            f.write(old_coord)
            f.writelines('\n')
    print("Old coordinates extraction complete.")
    return old_coords

#________________________________________________________________________________________________________________________________________________________

def get_sun_coordinates(ms, output_file):
    '''
    Extracts the coordinates of the Sun from a measurement set and writes them to a file.

    Parameters:
    ms (str): Path to the measurement set file.
    outfile (str): Path to the output file.

    Returns:
    None
    '''
    def format_coords(ra0,dec0):
        c = SkyCoord(ra0*u.deg,dec0*u.deg,frame='fk5')
        #hms = c.ra.to_string(u.hour, precision=1, pad=True)
        # Format Dec without decimal seconds
        #dms = c.dec.to_string(u.deg, precision=1, pad=True) 
        #return hms, dms
        hms = str(c.ra.to_string(u.hour, precision=1))
        dms = str(c.dec.to_string(u.deg, precision=1))
        #dms = str(c.dec)
        return hms,dms

    # MeerKAT
    obs_lat = -30.71323598930457
    obs_lon = 21.443001467965008
    loc = EarthLocation.from_geodetic(obs_lat, obs_lon)

    maintab = table(ms, ack=False)
    scans = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    lines = []

    print("Extracting Sun coordinates from {}...".format(ms))
    for scan in scans:
        subtab = maintab.query(query='SCAN_NUMBER==' + str(scan))
        t_scan = numpy.mean(subtab.getcol('TIME'))
        t = Time(t_scan / 86400.0, format='mjd')

        with solar_system_ephemeris.set('builtin'):
            sun = get_body('Sun', t, loc)
            sun_ra = sun.ra.value
            sun_dec = sun.dec.value
            sun_hms, sun_dms = format_coords(sun_ra, sun_dec)
            lines.append("{} {}".format(sun_hms, sun_dms))

    maintab.close()

    with open(output_file, 'wt') as f:
        for line in lines:
            f.write(line + '\n')

    print("Sun coordinates extracted and saved to {}.".format(output_file))


#________________________________________________________________________________________________________________________________________________________

def shift_coordinates(ms_list, coords, splitted_ms_dir, datacolumn='all'):
    '''
    A funtion that takes a list of scans and coordinates shift/rephase it for a specific colunm (CORRECTED_DATA) and iutput in the splitted_ms_dir directory 
    Parameters:
    ms_list (list): Path to the Measurement Sets.
    coords (File): Path to the coordinate file.
    splitted_ms_dir (Directory): Path to the scans directory
    datacolumn (str): Datacolumn to use (when not defined default is 'all').
    '''
    coordinates = []
    with open(coords, 'r') as file:
        for line in file:
            ra, dec = line.strip().split() # Assuming RA and Dec are separated by a space
            coordinates.append((ra,dec))

    #Sort the ms_list in numerical order
    sorted_ms_list = sorted(ms_list, key=lambda x: int(x.split('_scan_')[1].split('.')[0]))
    print(sorted_ms_list)
    chgcentre_path= '/home/samboco/solarKAT/Git_clone/wsclean/build/chgcentre'

    for ms, (ra, dec) in zip(sorted_ms_list, coordinates):
        #for ms in sorted_ms_list:
        command=[chgcentre_path, ms, ra, dec]
        try:
            subprocess.run(command, check=True)
            print("Successfully processed RA: {}, Dec:{} for MS: {}".format(ra, dec, ms))
        except subprocess.CalledProcessError as e:
            print("Error processing RA: {}, Dec: {} for MS: {}".format(ra, dec, ms))
            print("Error message: {}".format(e))

#________________________________________________________________________________________________________________________________________________________

def create_ds9_region_from_file(input_file, output_dir, ms):
    """
    Function to create a DS9 region file for each coordinate in the input file
    """
    #Open the MS table
    tab = table(ms, readonly=True)
    # Extract the scan numbers
    scan_numbers = list(np.unique(tab.getcol('SCAN_NUMBER')))
    # Close the MS table
    tab.close()

    with open(input_file, 'r') as f:
        for i, line in enumerate(f):
            ra, dec = map(str.strip, line.split())
            print('RA in hms format:', ra, 'DEC in dms format:', dec, output_dir)
            print('Processing coordinates for scan {}...'.format(scan_numbers[i]))

            ra_deg = hms2deg(ra)
            dec_deg = dms2deg(dec)
            print('RA in degrees:', ra_deg, 'DEC in degrees:', dec_deg)
            print('Creating DS9 region for scan {}...'.format(scan_numbers[i]))

            scan_number = scan_numbers[i]  # Use the scan number at the corresponding index
            print("Creating region for scan {} ............................................".format(scan_number))
            create_ds9_region(ra_deg, dec_deg, scan_number, output_dir)
    print("DS9 region creation completed successfully.")

def create_ds9_region(ra, dec, scan_number, output_dir):
    # Function to create a DS99 region file
    sun_region = f"""# Region file format: DS9 CARTA 3.0.0-beta.3 
global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle({ra}, {dec}, 1188.0000") # color=green
"""
#The UHF band images of the sun has a radius of 1062.9118 and the L bands 1188.0000"
    with open("{}/sun-region-{}.reg".format(output_dir, scan_number), 'w') as f:
        f.write(sun_region)

#_____________________________________________________________________________________________________________________________________________________

