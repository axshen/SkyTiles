#!/usr/bin/env python
import numpy
import csv
import os
import time
import argparse
import json
import sys



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("Generate tiles for a specfic SB.")
    parser.add_argument('-j', dest='json', help='The Json file.')
    args = parser.parse_args()
    config_json_file = args.json

    with open(config_json_file, "r") as read_file:
        tile_config = json.load(read_file)
    
    input_image = tile_config['input_image']

    input_fits_path  = input_image['path']
    input_image_name = input_image['name']
     
    imagefits = input_fits_path + input_image_name
    imExist = os.path.exists(imagefits)
    
    if not imExist:
        sys.exit('Input image or directory not found. Please check! Code aborted.')
    
    
    data_csv = tile_config['input_tile']
    
    csv_path = data_csv['path']
    csv_name = data_csv['csv_name']
    csv_file = csv_path + csv_name
    csvExist = os.path.exists(csv_file) 
    
    if not csvExist:
        sys.exit('>>>> CSV file or directory not found. Please check! Code aborted.')
    
    tile_template = tile_config['tile_template']
    
    template_path = tile_template['path']
    template_name = tile_template['name']
    template_file = template_path + template_name
    template_Exist = os.path.exists(template_file) 
    
    if not template_Exist:
        sys.exit('Template tile fits file or directory not found. Please check! Code aborted.')
    
    
    # make directories for the output tiles and read tile configuration.
    outputs = tile_config['outputs']
    
    naxis = outputs['naxis']
    parent_directory = outputs['path']
    output_prefix = outputs['tile_prefix']
    
    SB_ID   = csv_name.split('_')[-1].split('.')[0]   
    SB_directory = str(SB_ID) + '/' 

    tile_path = os.path.join(parent_directory, SB_directory)
    outdir_Exist = os.path.exists(tile_path)

    if not outdir_Exist:
        os.mkdir(tile_path)
        print("Output directory not found, so creating new ones named: '% s'." % tile_path)
    
    pixel_ID = []
    crpix1 = []
    crpix2 = []

    with open(csv_file, mode='r') as csv_info:
        csv_reader = csv.DictReader(csv_info)
    
        for row in csv_reader:
            pixel_ID.append(float(row['PIXELS']))
            crpix1.append(float(row['CRPIX_RA']))
            crpix2.append(float(row['CRPIX_DEC']))
       

    fitsheader = imhead(imagefits)
    axis = fitsheader['axisnames']
  
    #get template header
    template_header = imregrid(imagename=template_file, template='get')

    start_tiling = time.time()
    
    for i, (ra, dec) in enumerate(zip(crpix1, crpix2)):
        
        one_tile_start = time.time()
        
        #this is how I will update the dictionary
        template_header['csys']['direction0']['crpix'] = numpy.array([ra, dec])
    
        if len(axis) == 4:
            fourth_axis = axis[3]
            if fourth_axis == 'Frequency':
                number_of_frequency = fitsheader['shape'][3]
                template_header['shap'] = numpy.array([naxis, naxis, 1, number_of_frequency])
            
            third_axis = axis[2]
            if third_axis == 'Frequency':
                number_of_frequency = fitsheader['shape'][2]
                template_header['shap'] = numpy.array([naxis, naxis, number_of_frequency, 1])

        if len(axis) == 3:
            third_axis = axis[2]
            if third_axis == 'Frequency':
                number_of_frequency = fitsheader['shape'][2]
                template_header['shap'] = numpy.array([naxis, naxis, number_of_frequency]) 
            else:
                template_header['shap'] = numpy.array([naxis, naxis, 1])   

        outputname = tile_path + '%s-%s-%d.image'%(output_prefix, SB_ID, pixel_ID[i])
    
        # tiling, outputs tile fits in CASA image. 
        imregrid(imagename=imagefits, template=template_header, output=outputname, axes=[0, 1],  interpolation='cubic')
        #convert casa image to fits image
        one_tile_end = time.time()
        print('Tiling of pixel ID %d completed. Time elapsed %.3f seconds. '%(pixel_ID[i], (one_tile_end-one_tile_start)))
        
        print('Converting the casa image to fits image.')
        exportfits(imagename=outputname, fitsimage=outputname.split('.image')[0]+ '.fits')
        # delete all casa image files.
        print('Deleting the casa image. ')
        os.system('rm -rf %s'%outputname) 
        
    end_tiling = time.time()    
    print('Tiling for SB%s completed. Time elapsed is %.3f seconds.'%(SB_ID, (end_tiling - start_tiling)))
 

