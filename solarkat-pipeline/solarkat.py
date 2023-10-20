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
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body_barycentric, get_body, get_moon


#________________________________________________________________________________________________________________________________________________________

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

#________________________________________________________________________________________________________________________________________________________

def scan_numbers(ms, outfile):
    '''
    Extracts unique scan numbers from the Measurement Set
    ms (MS): Measurement Set
    output_file (str): Output file path to save the scan numbers
    '''
    from casacore.tables import table
    import numpy
    scans = []
    maintab = table(ms, ack=False)
    scan_no = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    for scan in scan_no:
        scans.append(str(scan))
    maintab.close()
    print(scans)
    
    # Save scan numbers to the output file
    with open(outfile, 'w') as file:
        file.write('\n'.join(scans))    
    return scans


def extract_scan_number(ms_scan):
    # Extract the scan number from the scan file name using regex
    scan_number = re.search(r"scan_(\d+)\.ms", ms_scan).group(1)
    return int(scan_number)
#________________________________________________________________________________________________________________________________________________________

def rename_model_data_column(ms_list, oldname, newname):
    '''
    Rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL' in a list of MS tables.
    Parameters:
    ms_list (str): Path to the measurement sets files.
    '''
    for ms_name in ms_list:
        # open the MS table in read-write mode
        ms = table(ms_name, readonly=False)
        # rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL'
        ms.renamecol(oldname, newname)
        # close the MS table
        ms.close()

#________________________________________________________________________________________________________________________________________________________

def perscan_timerange_intervals(interval, scan_list):
    '''
    A funtion that determines the number of intervals to devide the MS
    Parameters:
    interval (int): number of intervals on which to devide the measurement set
    Example: If interval=3 we the MS will be devided in three intervals (1-2, 2-3 and 3-4), which gives us three sub-scans of each scan
    scan_list (str): Path to the scan files
    ''' 
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


#________________________________________________________________________________________________________________________________________________________

def get_old_coords(ms_list, outfile):
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
        #print(old_coords)
        field_dir.close()

    with open(outfile,'wt') as f:
        for old_coord in old_coords:
            f.write(old_coord)
            f.writelines('\n')

    return old_coords

#________________________________________________________________________________________________________________________________________________________

def get_sun_coordinates(ms, outfile):
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

#________________________________________________________________________________________________________________________________________________________

def get_pertime_sun_coordinates(ms_list, outfile):
    '''
    Extracts the scan coordinates of the Sun from a list of  measurement sets splitted by time intervals and writes them to a file.

    Parameters:
    ms_list (str): Path to the measurement sets pertime files.
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
#________________________________________________________________________________________________________________________________________________________


def shift_coordinates(ms_list, coords, splitted_ms_dir):
    '''
    A funtion that takes a list of scans and coordinates shift/rephase it for a specific colunm (CORRECTED_DATA) and iutput in the splitted_ms_dir directory 
    Parameters:
    ms_list (str): Path to the measurement sets files.
    coords (File): Path to the coordinate files.
    splitted_ms_dir (Directory): Path to the directory
    datacolumn (str): Datacolumn to use
    '''
    with open(coords) as f:
        lines = f.readlines()
        num_lines = len(lines)

    #Sort the ms_list in numerical order
    sorted_ms_list = sorted(ms_list, key=lambda x: int(x.split('_scan_')[1].split('.')[0]))

    for i, ms in enumerate(sorted_ms_list):
        splitted_ms = os.path.join(splitted_ms_dir, os.path.basename(ms).replace(".ms", ".ms"))
        line = lines[i % num_lines].strip()
        #os.system(f"chgcentre {ms} {line} {splitted_ms} datacolumn={datacolumn}")
        print(line)
        os.system(f"chgcentre {ms} {line} {splitted_ms}")
        #subprocess.run(["chgcentre", ms, line, splitted_ms])
        print(f'Scan {ms} Done.')
    print("Phase center changed successfully.")


#________________________________________________________________________________________________________________________________________________________

def create_ds9_region_from_file(input_file, output_dir, ms):
    """
    Function to create a ds9 region file for each coordinate in the input file
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
            ra_deg = hms2deg(ra)
            dec_deg = dms2deg(dec)
            print('RA in degrees:', ra_deg, 'DEC in degrees:', dec_deg)
            scan_number = scan_numbers[i]  # Use the scan number at the corresponding index
            print(f"Creating region for scan {scan_number}"  '..................................................................................')
            create_ds9_region(ra_deg, dec_deg, scan_number, output_dir)

def create_ds9_region(ra, dec, scan_number, output_dir):
    # Function to create a ds9 region file
    sun_region = f'''# Region file format: DS9 CARTA 3.0.0-beta.3 
global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle({ra}, {dec}, 1188.0000") # color=green
'''
#The UHF band images of the sun has a radius of 1062.9118 and the L bands 1188.0000"
    with open(f'{output_dir}/sun_region_{scan_number}.reg', 'w') as f:
        f.write(sun_region)



# def create_ds9_region_from_file(input_file, output_dir, ms):
#     """
#     Function to create a ds9 region file for each coordinate in the input file
#     """
#     #Open the MS table
#     tab = table(ms, readonly=True)
#     # Extract the scan numbers
#     scan_numbers = list(np.unique(tab.getcol('SCAN_NUMBER')))
#     # Close the MS table
#     tab.close()

#     with open(input_file, 'r') as f:
#         for i, line in enumerate(f):
#             ra, dec = map(str.strip, line.split())
#             print('RA in hms format:', ra, 'DEC in dms format:', dec, output_dir)
#             ra_deg = hms2deg(ra)
#             dec_deg = dms2deg(dec)
#             print('RA in degrees:', ra_deg, 'DEC in degrees:', dec_deg)
#             scan_number = scan_numbers[i]  # Use the scan number at the corresponding index
#             print(f"Creating region for scan {scan_number}"  '..................................................................................')
#             if is_sun_present(ra_deg, dec_deg):  # Check if the Sun is present
#                 create_ds9_region(ra_deg, dec_deg, scan_number, output_dir)
#             else:
#                 print(f"Sun not detected in scan {scan_number}. Skipping region creation.")

# def create_ds9_region(ra, dec, scan_number, output_dir):
#     # Function to create a ds9 region file
#     sun_region = f'''# Region file format: DS9 CARTA 3.0.0-beta.3 
# global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
# fk5
# circle({ra}, {dec}, 1188.0000") # color=green
# '''
#     #The UHF band images of the sun has a radius of 1062.9118 and the L bands 1188.0000"
#     with open(f'{output_dir}/sun_region_{scan_number}.reg', 'w') as f:
#         f.write(sun_region)

# def is_sun_present(ra, dec):
#     # Function to check if the Sun is present in the image
#     # This is a placeholder function. You need to implement the actual Sun detection logic.
#     return True

#_____________________________________________________________________________________________________________________________________________________


def add_column_to_ms(ms, col_names, like_col):
    """
    Add columns to a single measurement file
    """
    success = False
    try:
        tb = table(ms, readonly=False)
    except Exception:
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
            tb.addcols(desc, dminfo)
            success = True
    
    tb.close()
    return success

#________________________________________________________________________________________________________________________________________________________

def add_columns_to_ms_list(ms_list, col_names, like_col):
    """
    Add columns to a list of measurement files
    """
    for ms_file in ms_list:
        try:
            tb = table(ms_file, readonly=False)
        except Exception:
            continue
        
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
                tb.addcols(desc, dminfo)
    
        tb.close()
#________________________________________________________________________________________________________________________________________________________


def rename_columns(ms_list):
    '''
    Rename columns in a list of MS tables.
    '''
    for ms_name in ms_list:
        '''
        Open the MS table in read-write mode
        '''
        ms = table(ms_name, readonly=False)
        '''
        Rename the 'MODEL_DATA' column to 'SUN_MODEL'
        '''
        if 'MODEL_DATA' in ms.colnames():
            ms.renamecol('MODEL_DATA', 'SUN_MODEL')
        # Rename the 'MODEL_DATA_ORIGINAL' column to 'MODEL_DATA'
        if 'MODEL_DATA_ORIGINAL' in ms.colnames():
            ms.renamecol('MODEL_DATA_ORIGINAL', 'MODEL_DATA')
        # close the MS table
        ms.close()

#________________________________________________________________________________________________________________________________________________________

def copy_columns(from_ms, to_ms, columns_to_copy):
    '''
This code defines a function copy_columns that takes in two parameters from_ms and to_ms which are the source and destination 
tables respectively. It also takes in a list of column names to copy from the source table to the destination table. 
The code then opens the two tables in a with block, reads the data from the source table for each column name in the list of 
columns to copy, and writes that data to the destination table.
''' 

    with table(from_ms, readonly=True) as src_table, table(to_ms, readonly=False) as dest_table:
        #columns_to_copy = ["MODEL_DATA", "MODEL_DATA_SUN"]
        for column_name in columns_to_copy:
            src_data = src_table.getcol(column_name)
            dest_table.putcol(column_name, src_data)

#________________________________________________________________________________________________________________________________________________________


#Copying Model data into the model-data-sun in the  ooriginal MS
def copy_model_data_to_model_data_sun(ms, ms_list, copycol, tocol):
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

    # Close the maintab
    maintab.close()

#________________________________________________________________________________________________________________________________________________________

#Copying a list of columns from the ms in ms-list into a list of columns in the original MS

def copy_columns_to_target(ms, ms_list, copycols, tocols):
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
            print("target rows", target_rows)

            # Get the columns' data from the source_subtab
            source_data = [source_subtab.getcol(copycol) for copycol in copycols]
            print("source data:", source_data)
            # Get the shapes of the source columns' data
            source_shapes = [data.shape for data in source_data]

            # Check if the target_subtab has the tocols columns
            for tocol, source_shape in zip(tocols, source_shapes):
                if tocol not in maintab.colnames():
                    # Add the tocol column to the target_subtab with the source_shape
                    maintab.addcols(
                        columns={tocol: {'datatype': 'complex', 'shape': source_shape}}
                    )

            # Get the columns' data from the target_subtab
            target_data = [target_rows.getcol(tocol) for tocol in tocols]
            print('Target data:', target_data)

            # Update the target columns' data with the source columns' data
            for target_col, source_col, data in zip(target_data, tocols, source_data):
                print(f"Copying column '{source_col}' to '{target_col}' or scan number {scan_number}...")
                target_col[:] = data

            # Put the updated columns' data back into the target_subtab
            for target_col in target_data:
                target_rows.putcol(tocol, target_col)

            # Flush the changes to disk
            maintab.flush()

            # Close the source_subtab
            source_subtab.close()

    # Close the maintab
    maintab.close()

#_________________________________________________________________________________________________________________________________________________________
##########################################################################################################################################################

import subprocess
from MSUtils.msutils import sumcols


# def subtract_model_data(ms_list):
#     #sorted_ms_list = sorted(ms_list)  # Sort the list of MS filenames
#     sorted_ms_list = sorted(ms_list, key=lambda x: int(x.split('_scan_')[1].split('.')[0]))

#     for ms in sorted_ms_list:
#         print("ms:", ms)
#         #residual = ms.query( "select CORRECTED_DATA  - MODEL_DATA as CORRECTED_DATA_SUN")

#         #print(residual)

#         #sumcols(ms, col1='MODEL_DATA_SUN', col2='CORRECTED_DATA', outcol= 'CORRECTED_DATA_SUN', subtract=True)
        
#         taql  "UPDATE {ms}  SET CORRECTED_DATA_SUN=CORRECTED_DATA-MODEL_DATA_SUN"

#subprocess.run(taql_command, shell=True)


def subtrac_sun(ms_list):
    sorted_ms_list = sorted(ms_list, key=lambda x: int(x.split('_scan_')[1].split('.')[0]))
    for ms in sorted_ms_list:
        print('ms:', ms)
        print('Subtracting CORRECTED_DATA_SUN from CORRECTED_DATA_SUN in', ms)
        os.system(f"/usr/bin/taql update {ms} set CORRECTED_DATA_SUN=CORRECTED_DATA-MODEL_DATA_SUN")
    
