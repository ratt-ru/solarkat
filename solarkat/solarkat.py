#!/usr/bin/env python
"""
Helper functions called by the SolarKAT stimela recipe (solarkat.yaml) through
the cabs defined in solarkat-cabs.yaml.
 
Coordinate files written here have one line per scan:
 
    <scan_number> <ra> <dec>
 
with RA as 00h00m00.0s and Dec as -00d00m00.0s, which is the format the
chgcentre cab (cult-cargo) accepts. Keeping the scan number on each line means
that a coordinate is always matched to its scan by number, never by position.
"""
import os
import re
 
import numpy as np
from astropy import units as u
from astropy.coordinates import (GCRS, ICRS, Angle, EarthLocation, SkyCoord,
                                 get_body, solar_system_ephemeris)
from astropy.time import Time
from casacore.tables import table
 
# MeerKAT array reference position
MEERKAT_LAT_DEG = -30.71323598930457
MEERKAT_LON_DEG = 21.443001467965008
 
# Default radius of the DS9 region (and therefore of the clean mask) drawn around
# the Sun, in arcsec. These are the values used so far; override them with the
# recipe's sun_radius_arcsec setting.
DEFAULT_SUN_RADIUS_ARCSEC = {
    "L": 1188.0,
    "UHF": 1062.9118,
}
 
# ______________________________________________________________________________
# Auxiliary functions
 
 
def hms2deg(hms):
    '''Convert an hms string to degrees.'''
    return Angle(hms, unit='hour').degree
 
 
def dms2deg(dms):
    '''Convert a dms string to degrees.'''
    return Angle(dms, unit='degree').degree
 
 
def extract_scan_number(ms_scan):
    '''Return the scan number of a per-scan MS named like <name>_scan_<N>.ms.'''
    match = re.search(r"_scan_(\d+)\.ms$", os.path.basename(os.path.normpath(ms_scan)))
    if match is None:
        raise ValueError("Cannot find a scan number in MS name {}".format(ms_scan))
    return int(match.group(1))
 
 
def format_ra_dec(ra_deg, dec_deg, precision=3):
    '''Format RA/Dec (degrees) in the 00h00m00.0s / -00d00m00.0s form chgcentre accepts.'''
    c = SkyCoord(ra_deg * u.deg, dec_deg * u.deg, frame='icrs')
    ra = c.ra.to_string(unit=u.hourangle, sep='hms', pad=True, precision=precision)
    dec = c.dec.to_string(unit=u.degree, sep='dms', pad=True, precision=precision)
    return ra, dec
 
 
def write_coords_file(output_file, rows):
    '''Write (scan, ra, dec) rows to a coordinate file.'''
    with open(output_file, 'wt') as f:
        for scan, ra, dec in rows:
            f.write("{} {} {}\n".format(scan, ra, dec))
 
 
def read_coords_file(coords_file):
    '''Read a coordinate file into an ordered dict {scan(str): (ra, dec)}.'''
    coords = {}
    with open(coords_file) as f:
        for line in f:
            if not line.strip():
                continue
            scan, ra, dec = line.split()
            coords[scan] = (ra, dec)
    return coords
 
 
# ______________________________________________________________________________
 
 
def extract_and_save_scan_numbers(ms, output_file):
    '''
    Extracts unique scan numbers from the Measurement Set and saves them, one per line.
 
    ms (str): Measurement Set path
    output_file (str): Output file path to save the scan numbers
    Returns:
    List[str]: List of unique scan numbers as strings
    '''
    print("Extracting scan numbers from {}...".format(ms))
    with table(ms, ack=False) as maintab:
        scans = [str(scan) for scan in np.unique(maintab.getcol('SCAN_NUMBER'))]
 
    print("Unique scan numbers extracted:", scans)
    with open(output_file, 'w') as f:
        f.write('\n'.join(scans) + '\n')
    print("Scan numbers saved to {}".format(output_file))
    return scans
 
 
# ______________________________________________________________________________
 
 
def get_old_coords(ms, output_file):
    '''
    Writes the original phase centre of every scan to a coordinate file.
 
    The phase centres are read from the FIELD table of the original MS (not from
    the per-scan MSs), so the result stays correct if the recipe is re-run after
    the per-scan MSs have already been shifted to the Sun.
 
    ms (str): Path to the original Measurement Set.
    output_file (str): Path to the output coordinate file.
    '''
    with table(ms, ack=False) as maintab:
        scan_col = maintab.getcol('SCAN_NUMBER')
        field_col = maintab.getcol('FIELD_ID')
    with table("{}::FIELD".format(ms), ack=False) as field_tab:
        phase_dir = field_tab.getcol("PHASE_DIR")  # radians, shape (nfield, npoly, 2)
 
    rows = []
    for scan in np.unique(scan_col):
        fields = np.unique(field_col[scan_col == scan])
        if len(fields) != 1:
            raise RuntimeError("Scan {} contains {} fields; expected one.".format(scan, len(fields)))
        ra_rad, dec_rad = phase_dir[fields[0], 0]
        ra, dec = format_ra_dec(np.degrees(ra_rad), np.degrees(dec_rad))
        rows.append((int(scan), ra, dec))
        print("Scan {}: original phase centre {} {}".format(scan, ra, dec))
 
    write_coords_file(output_file, rows)
    print("Original phase centres saved to {}".format(output_file))
 
 
# ______________________________________________________________________________
 
 
def sun_direction(t, location):
    '''
    Direction of the Sun to use as an ICRS/J2000 phase centre at time t.
 
    get_body() returns the apparent topocentric position in GCRS. Casacore treats a
    J2000 phase centre like a distant source and adds aberration itself, so we
    convert the apparent direction back to ICRS without a distance; this removes
    the ~20 arcsec of annual aberration that would otherwise be applied twice.
    '''
    with solar_system_ephemeris.set('builtin'):
        sun = get_body('Sun', t, location)
    apparent = SkyCoord(sun.ra, sun.dec,
                        frame=GCRS(obstime=t, obsgeoloc=sun.obsgeoloc, obsgeovel=sun.obsgeovel))
    icrs = apparent.transform_to(ICRS())
    return icrs.ra.deg, icrs.dec.deg
 
 
def get_sun_coordinates(ms, output_file):
    '''
    Computes the position of the Sun at the mean time of every scan and writes
    them to a coordinate file.
 
    ms (str): Path to the measurement set file.
    output_file (str): Path to the output coordinate file.
    '''
    # from_geodetic takes (lon, lat)
    loc = EarthLocation.from_geodetic(lon=MEERKAT_LON_DEG * u.deg, lat=MEERKAT_LAT_DEG * u.deg)
 
    with table(ms, ack=False) as maintab:
        scan_col = maintab.getcol('SCAN_NUMBER')
        time_col = maintab.getcol('TIME')  # MJD seconds, UTC
 
    print("Extracting Sun coordinates from {}...".format(ms))
    rows = []
    for scan in np.unique(scan_col):
        t = Time(np.mean(time_col[scan_col == scan]) / 86400.0, format='mjd', scale='utc')
        ra, dec = format_ra_dec(*sun_direction(t, loc), precision=2)
        rows.append((int(scan), ra, dec))
        print("Scan {} ({}): Sun at {} {}".format(scan, t.isot, ra, dec))
 
    write_coords_file(output_file, rows)
    print("Sun coordinates saved to {}".format(output_file))
 
 
# ______________________________________________________________________________
 
 
def create_ds9_region_from_file(input_file, output_dir, band='L', radius_arcsec=None):
    """
    Creates a DS9 region file (sun-region-<scan>.reg) around the Sun for every scan
    in the coordinate file.
 
    band (str): observing band, used to pick a default radius.
    radius_arcsec (float): radius of the region; overrides the band default.
    """
    if radius_arcsec is None:
        if band not in DEFAULT_SUN_RADIUS_ARCSEC:
            raise ValueError("No default Sun radius for band '{}'; set sun_radius_arcsec.".format(band))
        radius_arcsec = DEFAULT_SUN_RADIUS_ARCSEC[band]
 
    os.makedirs(output_dir, exist_ok=True)
    for scan, (ra, dec) in read_coords_file(input_file).items():
        ra_deg, dec_deg = hms2deg(ra), dms2deg(dec)
        print("Scan {}: region at RA {} deg, Dec {} deg, radius {} arcsec".format(
            scan, ra_deg, dec_deg, radius_arcsec))
        create_ds9_region(ra_deg, dec_deg, scan, output_dir, radius_arcsec)
    print("DS9 region creation completed successfully.")
 
 
def create_ds9_region(ra, dec, scan_number, output_dir, radius_arcsec):
    '''Write one circular DS9 region file.'''
    sun_region = f"""# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle({ra}, {dec}, {radius_arcsec}") # color=green
"""
    with open(os.path.join(output_dir, "sun-region-{}.reg".format(scan_number)), 'w') as f:
        f.write(sun_region)
 
 
# ______________________________________________________________________________
 
 
def copy_scan_column_to_ms(ms, ms_list, fromcol, tocol):
    '''
    Copies column `fromcol` of every per-scan MS into column `tocol` of the
    original MS, on the rows of the matching scan.
 
    The per-scan MSs come from CASA split with no averaging, so each holds the
    same rows, in the same order, as its scan in the original MS. This is checked
    on TIME, ANTENNA1, ANTENNA2 and DATA_DESC_ID before anything is written.
    '''
    scan_ms = {extract_scan_number(m): m for m in ms_list}
    check_cols = ['TIME', 'ANTENNA1', 'ANTENNA2', 'DATA_DESC_ID']
 
    with table(ms, readonly=False, ack=False) as maintab:
        if tocol not in maintab.colnames():
            raise RuntimeError("Column {} does not exist in {}; add it first.".format(tocol, ms))
        scans = np.unique(maintab.getcol('SCAN_NUMBER'))
        missing = [int(s) for s in scans if int(s) not in scan_ms]
        if missing:
            raise RuntimeError("No per-scan MS found for scan(s) {}".format(missing))
 
        for scan in scans:
            with table(scan_ms[int(scan)], ack=False) as src, \
                    maintab.query('SCAN_NUMBER=={}'.format(int(scan))) as dst:
                if src.nrows() != dst.nrows():
                    raise RuntimeError("Scan {}: {} rows in {} but {} in {}".format(
                        scan, src.nrows(), scan_ms[int(scan)], dst.nrows(), ms))
                for col in check_cols:
                    if not np.array_equal(src.getcol(col), dst.getcol(col)):
                        raise RuntimeError("Scan {}: column {} differs between {} and {}".format(
                            scan, col, scan_ms[int(scan)], ms))
                dst.putcol(tocol, src.getcol(fromcol))
            print("Copied {} -> {} for scan {}".format(fromcol, tocol, scan))
        maintab.flush()
    print("Copying {} to {} in {} completed.".format(fromcol, tocol, ms))
 
