#!/usr/bin/env python
from astropy.io import fits
import json
import pandas as pd
import os
import sys
import argparse



"""

    This script generates tiles for a specified field.
    Inputs:
	1. Json file containing all the information to the needed input files.   
        2. CSV file for the corresponding field. This file is generated using PyMapSkyTiles.py.
           It contains the CRPIXs needed for tilling as well as the tile IDs. The tile name 
           is 'filename-SBID.csv' where SBID specifies the ID of an SB e.g. '2156-54'.
           NB: Parse via the json file under 'input_tile_config'.
        3. Image MFS/cube image. Parse via the json file using 'input_image'.

    Output:
        1. Headers (one for each tile) '.hdr'.
        2. Tiles.   
        NB: Specific the prefix of the output tile images and headers using 'output_tile'. 
     
Author: Lerato

"""


def Tile_Header_Cube(ra_pixel, dec_pixel, cdelt, naxis, crval_ra, crval_dec):

    """ Headers for the corresponding tiles. Used during the reprojection.

        ra_pixel/dec_pixel : CRPIX1/CRPIX2.
        cdelt              : Pixel size of the tiles.
        naxis              : Tile size in pixels.
        crval_ra/crval_dec : CRVAL1/CRVAL2 in degs. Written as part of the comments for future reference.

    """

    HPX_2D_header  =  "SIMPLE  =          T / file does conform to FITS standard \n"
    HPX_2D_header +=  "BITPIX  =          -32 / number of bits per data pixel \n"
    HPX_2D_header +=  "NAXIS   =           3 / number of data axes \n"
    HPX_2D_header +=  "NAXIS1  =           %d / length of data axis 1 \n" %naxis
    HPX_2D_header +=  "NAXIS2  =           %d / length of data axis 2 \n" %naxis
    HPX_2D_header +=  "NAXIS3  =           %d / length of data axis 2 \n" %no_of_frequencies
    HPX_2D_header +=  "BUNIT   = 'Jy/beam'   / Brightness (pixel) unit \n"
    HPX_2D_header +=  "EQUINOX =   2.000000000000E+03 \n"
    HPX_2D_header +=  "RADESYS = 'FK5' \n"    
    HPX_2D_header +=  "PC1_1   =           0.70710677 / Transformation matrix element \n"
    HPX_2D_header +=  "PC1_2   =           0.70710677 / Transformation matrix element \n"
    HPX_2D_header +=  "PC2_1   =           -0.70710677 / Transformation matrix element \n"
    HPX_2D_header +=  "PC2_2   =           0.70710677 / Transformation matrix element \n"
    HPX_2D_header +=  "CRPIX1  =           %r / Coordinate reference pixel \n" %ra_pixel
    HPX_2D_header +=  "CRPIX2  =           %r / Coordinate reference pixel \n" %dec_pixel
    HPX_2D_header +=  "CDELT1  =          -%r / [deg] Coordinate increment \n" %cdelt
    HPX_2D_header +=  "CDELT2  =           %r / [deg] Coordinate increment \n" %cdelt
    HPX_2D_header +=  "CTYPE1  = 'RA---HPX'  / Right ascension in an HPX projection \n" 
    HPX_2D_header +=  "CTYPE2  = 'DEC--HPX'  / Declination in an HPX projection \n"
    HPX_2D_header +=  "CRVAL1  =             0. / [deg] Right ascension at the reference point \n" 
    HPX_2D_header +=  "CRVAL2  =             0. / [deg] Declination at the reference point \n" 
    HPX_2D_header +=  "CUNIT1  = 'deg     ' \n"
    HPX_2D_header +=  "CUNIT2  = 'deg     ' \n"
    HPX_2D_header +=  "CTYPE3  =  %r \n"%freq_name
    HPX_2D_header +=  "CRVAL3  =   %s \n"%("{:.9E}".format(freq_cval))                                                  
    HPX_2D_header +=  "CDELT3  =   %s \n"%("{:.9E}".format(freq_cdelt))
    HPX_2D_header +=  "CRPIX3  =   %s \n"%("{:.9E}".format(freq_cpix))
    HPX_2D_header +=  "CUNIT3  = 'Hz    ' \n"
    HPX_2D_header +=  "CTYPE4  = %r \n"%stokes_name
    HPX_2D_header +=  "CRVAL4  =   %s \n" %("{:.9E}".format(stokes_cval))
    HPX_2D_header +=  "CDELT4  =   %s \n" %("{:.9E}".format(stokes_cdelt))
    HPX_2D_header +=  "CRPIX4  =   %s \n" %("{:.9E}".format(stokes_cpix))
    HPX_2D_header +=  "PV2_1   =             4 / HPX H parameter (longitude) \n" 
    HPX_2D_header +=  "PV2_2   =             3 / HPX K parameter  (latitude) \n"
    HPX_2D_header +=  "BMAJ    = %r\n"%bmaj
    HPX_2D_header +=  "BMIN    = %r\n"%bmin
    HPX_2D_header +=  "BPA     = %r\n"%bpa
    HPX_2D_header +=  "COMMENT CRVAL1 = %s degrees \n"%("{:.9E}".format(crval_ra))
    HPX_2D_header +=  "COMMENT CRVAL2 = %s degrees \n"%("{:.9E}".format(crval_dec))
    HPX_2D_header +=  "END" 

    return HPX_2D_header



def edit_data_header(fits_data, header):

    """ Edit the header and reduce the data for situations where the 
    input image has four axis (naxis)."""

    ctype3 = header['CTYPE3']
    ctype4 = header['CTYPE4']


    new_header = header.copy()
    new_header['naxis'] = 3

    if ctype4 == 'FREQ':
        new_header['ctype3'] = ctype4
        new_header['crval3'] = header['crval4']
        new_header['cunit3'] = header['cunit4']
        new_header['crpix3'] = header['crpix4']
        new_header['naxis3'] = header['naxis4']

        new_header['ctype4'] = ctype3
        new_header['crval4'] = header['crval3']
        new_header['cunit4'] = header['cunit3']
        new_header['crpix4'] = header['crpix3']
        fits_data = fits_data[:, 0, ...]
        new_header.remove('naxis4')

    if ctype3 == 'FREQ':
        fits_data = fits_data[0, ...]
        new_header.remove('naxis4')

    return fits_data, new_header


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("Generate tiles for a specfic SB.")
    parser.add_argument('-j', dest='json', help='The Json file.')
    args = parser.parse_args()
    config_json_file = args.json

    with open(config_json_file, "r") as read_file:
        tile_config = json.load(read_file)

    
    input_image = tile_config['input_image']

    fits_path  = input_image['path']
    image_name = input_image['name']
    image_file = fits_path + image_name

    # open the fits file.
    image  =  fits.open(image_file)
    fits_header = image[0].header
    fits_naxis = fits_header['naxis']

    global freq_name, freq_cval, freq_cdelt, freq_cpix, freq_unit, stokes_name, stokes_cval, stokes_cpix,\
        stokes_cdelt, no_of_frequencies, bmaj, bmin, bpa
    

    bmaj = fits_header['BMAJ']
    bmin = fits_header['BMIN']
    bpa  = fits_header['BPA']
   
    ctype3 = fits_header['ctype3']
    crval3 = fits_header['crval3']
    cdelt3 = fits_header['cdelt3']
    crpix3 = fits_header['crpix3']
    cunit3 = fits_header['cunit3']

    ctype4 = fits_header['ctype4']
    crval4 = fits_header['crval4']
    cdelt4 = fits_header['cdelt4']
    crpix4 = fits_header['crpix4']
    cunit4 = fits_header['cunit4']

    if ctype3 == 'STOKES':
        stokes_cval = crval3
        stokes_name = ctype3
        stokes_cpix = crpix3
        stokes_cdelt = cdelt3
    
        no_of_frequencies = fits_header['naxis4']
        freq_cval  = crval4
        freq_name  = ctype4
        freq_cpix  = crpix4
        freq_cdelt = cdelt4
        freq_unit = cunit4
    
    
    if ctype4 == 'STOKES':
        stokes_cval  = crval4
        stokes_name  = ctype4
        stokes_cpix  = crpix4
        stokes_cdelt = cdelt4
    
        no_of_frequencies = fits_header['naxis3']
        freq_cval  = crval3
        freq_name  = ctype3
        freq_cpix  = crpix3
        freq_cdelt = cdelt3
        freq_unit  = cunit3

    output_tile = tile_config['output_tile']
    naxis = output_tile['naxis']
    cdelt = output_tile['cdelt']
    outdir_tile = output_tile['path']
    outprefix_fitstiles = output_tile['fits_prefix']

    if stokes_cval == 1.0:
        outprefix_fitstiles = outprefix_fitstiles + '-I'
        stokes = 'I'
        
    if stokes_cval == 2.0: 
        outprefix_fitstiles = outprefix_fitstiles + '-Q'
        stokes = 'Q'
        
    if stokes_cval == 3.0:
        outprefix_fitstiles = outprefix_fitstiles + '-U'
        stokes = 'U'
        
    if stokes_cval == 4.0: 
        outprefix_fitstiles = outprefix_fitstiles + '-V'
        stokes = 'V'


    data_excel = tile_config['input_tile_config']
    excel_files_path = data_excel['path']
    config_SB_tile = data_excel['SB_tile']
    SB_ID   = config_SB_tile.split('_')[-1].split('.')[0]

    if excel_files_path is None:
        try: 
            filename = config_SB_tile
            df =  pd.read_csv(filename)
            print(df) 
        except:
            sys.exit('Excel files not found and no path/directory is provided.')
    else:
        filename = excel_files_path + config_SB_tile
        try: 
            df =  pd.read_csv(filename)
            print(df) 
        except:
            sys.exit('No excel file inside the provided path/directory.')
            
    
    healpix_pixels = [p for p in df['PIXELS']]
    CRVAL_RA  = [p for p in df['CRVAL_RA [deg]']]
    CRVAL_DEC = [p for p in df['CRVAL_DEC [deg]']]
    CRPIX_RA  = [p for p in df['CRPIX_RA']]
    CRPIX_DEC = [p for p in df['CRPIX_DEC']]
    

    if no_of_frequencies == 1:
        print('>>> Image provided has 1 channel, so it will be treated as an MFS, projected using mProject.')
    if no_of_frequencies > 1:
        print('>>> You provided an image CUBE. It will be processed using mprojectCube')
    
    # fix situations where the input image has 4 axis. Montage is unable 
    # to deal with such situations.
    if fits_naxis > 3: 
            
        print('>>> The input image data has %d number of axis. This will cause problems with Montage. '
                  'We going to attempt to reduce the image to 3D: (ra, dec,frequency).'%fits_naxis)
            
        newfits_data, newfits_hdr = edit_data_header(image[0].data, fits_header)
        tempfits_file = fits_path + 'temp-%s-%s.fits'%(SB_ID, stokes)
        fits.writeto(tempfits_file, newfits_data, newfits_hdr, overwrite=True)
            
        # set the image file to be the temporary file
        image_file = tempfits_file
    
    for i, (ra, dec) in enumerate(zip(CRPIX_RA, CRPIX_DEC)):
    
    
        header =   Tile_Header_Cube(ra, dec, cdelt, naxis, CRVAL_RA[i], CRVAL_DEC[i])
        header_output_file = outdir_tile + '%s-%s-header-%d.hdr'%(outprefix_fitstiles,SB_ID, healpix_pixels[i])
    
        header_file = open(header_output_file, 'w')
        header_file.write(str(header))
        header_file.close()
        
           
        if no_of_frequencies == 1:
            
            output_file = outdir_tile + '%s-%s-tile-%d.fits'%(outprefix_fitstiles, SB_ID, healpix_pixels[i])
            execute_string = 'mProject %s %s %s -f'%(image_file, output_file, header_file.name)
            os.system(execute_string)
            print(execute_string)
              
        if no_of_frequencies > 1:
            output_file = outdir_tile + '%s-%s-tile-%d.fits'%(outprefix_fitstiles, SB_ID, healpix_pixels[i])
            execute_string = 'mProjectCube %s %s %s -f'%(image_file, output_file, header_file.name)
            os.system(execute_string)
            print(execute_string)
                   
        
    if fits_naxis > 3:
        os.system('rm -rf %s'%tempfits_file)
