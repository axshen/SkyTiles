#!/usr/bin/env python3

import os
import sys
import numpy
import csv
import time
import argparse
import json
from casatasks import imhead, imregrid, exportfits


def main(argv):
    # Load configuration
    parser = argparse.ArgumentParser("Generate tiles for a specfic SB.")
    parser.add_argument('-j', dest='json', help='The JSON configuration file.')
    args = parser.parse_args(argv)
    with open(args.json, "r") as read_file:
        tile_config = json.load(read_file)

    # Read image, sky tiling map and tile template
    image = tile_config['image']
    if not os.path.exists(image):
        raise Exception(f'Input image {image} not found.')

    tiling_map = tile_config['tiling_map']
    if not os.path.exists(tiling_map):
        raise Exception(f'Input sky tiling map {tiling_map} not found.')

    tile_template = tile_config['tile_template']
    if not os.path.exists(tile_template):
        raise Exception(f'Input tile template {tile_template} not found.')

    # Outputs
    outputs = tile_config['outputs']
    naxis = outputs['naxis']
    output_prefix = outputs['tile_prefix']
    work_dir = outputs['path']
    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    sbid = tiling_map.split('/')[-1].split('_')[-1].split('.')[0]
    write_dir = os.path.join(work_dir, sbid)
    if not os.path.exists(write_dir):
        os.mkdir(write_dir)
        print(f"Output directory not found. Creating new directory: {write_dir}")

    # Read tiling map
    crpix1 = []
    pixel_ID = []
    crpix2 = []

    with open(tiling_map, mode='r') as f:
        csv_reader = csv.DictReader(f)
        for row in csv_reader:
            pixel_ID.append(float(row['PIXELS']))
            crpix1.append(float(row['CRPIX_RA']))
            crpix2.append(float(row['CRPIX_DEC']))

    fitsheader = imhead(image)
    axis = fitsheader['axisnames']

    # Read tile template header
    template_header = imregrid(imagename=tile_template, template='get')

    # Starting the tiling.
    start_tiling = time.time()
    for i, (ra, dec) in enumerate(zip(crpix1, crpix2)):
        one_tile_start = time.time()

        # this is how I will update the dictionary
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

        outputname = work_dir + '%s-%s-%d.image'%(output_prefix, sbid, pixel_ID[i])

        # tiling, outputs tile fits in CASA image.
        imregrid(imagename=image, template=template_header, output=outputname, axes=[0, 1],  interpolation='cubic')
        #convert casa image to fits image
        one_tile_end = time.time()
        print('Tiling of pixel ID %d completed. Time elapsed %.3f seconds. '%(pixel_ID[i], (one_tile_end - one_tile_start)))

        print('Converting the casa image to fits image.')
        exportfits(imagename=outputname, fitsimage=outputname.split('.image')[0]+ '.fits')
        # delete all casa image files.
        print('Deleting the casa image. ')
        os.system('rm -rf %s' % outputname)

    end_tiling = time.time()
    print('Tiling for SB%s completed. Time elapsed is %.3f seconds.'%(sbid, (end_tiling - start_tiling)))



if __name__ == '__main__':
    main(sys.argv[1:])
