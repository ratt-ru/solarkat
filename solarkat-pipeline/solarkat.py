#!/usr/bin/env python
import re
import numpy
import sys, os
import subprocess
import numpy as np
from astropy.time import Time
from pyrap.tables import table
from astropy import units as u
from casacore.tables import table
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
    print(f"Extracting scan numbers from {ms}...")
    
    scans = []
    
    with table(ms, ack=False) as maintab:
        scan_numbers = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
        scans = [str(scan) for scan in scan_numbers]

    print("Unique scan numbers extracted:")
    print(scans)

    # Save scan numbers to the output file
    with open(output_file, 'w') as file:
        file.write('\n'.join(scans))

    print(f"Scan numbers saved to {output_file}")
    
    return scans


#________________________________________________________________________________________________________________________________________________________

def rename_model_data_column(ms, oldname, newname):
    
    '''
    Rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL' in a Measurement Set (MS) table.
    Parameters:
    ms (str): Path to the Measurement Set file.
    oldname (str): Name of the column to be renamed.
    newname (str): New name for the column.
    '''
    print(f"Processing {ms}...")
    # open the MS table in read-write mode
    ms = table(ms, readonly=False)
    # rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL'
    ms.renamecol(oldname, newname)
    print(f"Column '{oldname}' renamed to '{newname}' in {ms}.")
    # close the MS table
    ms.close()

##########################################################################################################################################################

def rename_columns(ms_list, oldname, newname):
    '''
    Rename columns in a list of MS tables.
    '''
    for ms_name in ms_list:

        print(f"Renaming column '{oldname}' to '{newname}' in {ms_name}...")

        '''
        Open the MS table in read-write mode
        '''
        ms = table(ms_name, readonly=False)

         # Rename the 'oldname' column to 'newname'
        if 'oldname' in ms.colnames():
            ms.renamecol('oldname', 'newname')
        # close the MS table
        ms.close()
        print(f"Renaming completed for {ms_name}.")
    print("Rename column completed successfully.")

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

        print(f"Processing {ms}...")
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
        field_dir.close()
        print(f"Processing of {ms} complete.")

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

    print(f"Extracting Sun coordinates from {ms}...")
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

    with open(output_file, 'wt') as f:
        for line in lines:
            f.write(line + '\n')

    print(f"Sun coordinates extracted and saved to {output_file}.")


#________________________________________________________________________________________________________________________________________________________


def shift_coordinates(ms_list, coords, splitted_ms_dir):
    '''
    A funtion that takes a list of scans and coordinates shift/rephase it for a specific colunm (CORRECTED_DATA) and iutput in the splitted_ms_dir directory 
    Parameters:
    ms_list (list): Path to the Measurement Sets.
    coords (File): Path to the coordinate file.
    splitted_ms_dir (Directory): Path to the scans directory
    datacolumn (str): Datacolumn to use (when not defined default is 'all').
    '''
    with open(coords) as f:
        lines = f.readlines()
        num_lines = len(lines)

    #Sort the ms_list in numerical order
    sorted_ms_list = sorted(ms_list, key=lambda x: int(x.split('_scan_')[1].split('.')[0]))

    for i, ms in enumerate(sorted_ms_list):
        splitted_ms = os.path.join(splitted_ms_dir, os.path.basename(ms).replace(".ms", ".ms"))
        line = lines[i % num_lines].strip()
        print(f"Changing phase centre coordinates for {ms}...")
        os.system(f"chgcentre {ms} {line} {splitted_ms}")
        print(f'Scan {ms} Done.')
    print("Phase center changed successfully.")


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
            print(f'Processing coordinates for scan {scan_numbers[i]}...')

            ra_deg = hms2deg(ra)
            dec_deg = dms2deg(dec)
            print('RA in degrees:', ra_deg, 'DEC in degrees:', dec_deg)
            print(f'Creating DS9 region for scan {scan_numbers[i]}...')

            scan_number = scan_numbers[i]  # Use the scan number at the corresponding index
            print(f"Creating region for scan {scan_number}"  '..................................................................................')
            create_ds9_region(ra_deg, dec_deg, scan_number, output_dir)
    print("DS9 region creation completed successfully.")

def create_ds9_region(ra, dec, scan_number, output_dir):
    # Function to create a DS99 region file
    sun_region = f'''# Region file format: DS9 CARTA 3.0.0-beta.3 
global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle({ra}, {dec}, 1188.0000") # color=green
'''
#The UHF band images of the sun has a radius of 1062.9118 and the L bands 1188.0000"
    with open(f'{output_dir}/sun_region_{scan_number}.reg', 'w') as f:
        f.write(sun_region)

#_____________________________________________________________________________________________________________________________________________________


def add_column_to_ms(ms, col_names, like_col):
    """
    Add columns to a single measurement file
    """
    success = False
    try:
        print(f"Opening {ms}...")
        tb = table(ms, readonly=False)
    except Exception as e:
        print(f"Error: {e}")
        return success
    
    for col_name in col_names:
        if col_name not in tb.colnames():
            """
            Get column description from column 'like_col'
            """
            desc = tb.getcoldesc(like_col)
            desc[str('name')] = str(col_name)
            desc[str('comment')] = str(desc['comment'].replace(" ", "_"))
            dminfo = tb.getdminfo(like_col)
            dminfo[str("NAME")] = "{}-{}".format(dminfo["NAME"], col_name) 
            print(f"Adding column '{col_name}' to {ms}...")       
            tb.addcols(desc, dminfo)
            success = True
    print(f"Columns added to {ms} successfully.")
    tb.close()
    return success


#________________________________________________________________________________________________________________________________________________________


#Copying Model data into the model-data-sun in the  ooriginal MS

def copy_model_data_to_model_data_sun(ms, ms_list, copycol, tocol):
    print(f"Copying {copycol} to {tocol} in {ms}...")

    # Open the main MS table as maintab
    maintab = table(ms, readonly=False)

    # Get the unique scan numbers from the SCAN_NUMBER column in maintab
    scans = list(np.unique(maintab.getcol('SCAN_NUMBER')))

    for scan_number in scans:

        # Find the corresponding MS scan file for the current scan number
        ms_scan = next((ms_scan for ms_scan in ms_list if extract_scan_number(ms_scan) == scan_number), None)

        if ms_scan:
            # Open the MS scan as the source_subtab
            source_subtab = table(ms_scan, readonly=True)

            # Query the target_subtab to get the rows for the current scan number
            target_rows = maintab.query(query='SCAN_NUMBER==' + str(scan_number))

            # Get the MODEL_DATA from the source_subtab
            source_model_data = source_subtab.getcol(copycol)

            # Get the shape of the source_model_data
            source_shape = source_model_data.shape

            # Check if the target_subtab has the tocol column
            if tocol not in target_rows.colnames():
                # Add the tocol column to the target_subtab with the source_shape
                target_rows.addcols(
                    columns={tocol: {'datatype': 'complex', 'shape': source_shape}}
                )

            # Get the MODEL_DATA_SUN from the target_subtab
            target_model_data_sun = target_rows.getcol(tocol)

            # Update the target_model_data_sun with the source_model_data
            target_model_data_sun[:] = source_model_data

            # Put the updated MODEL_DATA_SUN back into the target_subtab
            target_rows.putcol(tocol, target_model_data_sun)

            # Flush the changes to disk
            maintab.flush()

            # Close the source_subtab
            source_subtab.close()
            print(f"Copy completed for scan {scan_number}.")

    # Close the maintab
    maintab.close()

    print(f"Copying {copycol} to {tocol} in {ms} completed successfully.")


#________________________________________________________________________________________________________________________________________________________

