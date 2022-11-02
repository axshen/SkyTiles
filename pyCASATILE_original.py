import numpy
import csv
import os
from casatasks import imhead, imregrid, importfits, exportfits


tile_template = 'data/tile-template.fits'
imagefits = 'data/image.restored.i.SB10040-Single.fits'
csv_file = 'PoSSUM_2156-54.csv'
outdir = 'tiles/MFS/'


pixel_ID = []
crpix1 = []
crpix2 = []

with open(csv_file, mode='r') as csv_info:
    csv_reader = csv.DictReader(csv_info)

    for row in csv_reader:
        pixel_ID.append(float(row['PIXELS']))
        crpix1.append(float(row['CRPIX_RA']))
        crpix2.append(float(row['CRPIX_DEC']))


imagecasa = imagefits.split('.fits')[0] + '.image'

fitsheader = imhead(imagefits)
axis = fitsheader['axisnames']

#convert the fits image to casa image
importfits(imagename=imagecasa, fitsimage=imagefits, overwrite=True)

#get template header
template_header = imregrid(imagename=tile_template, template='get')


SB_ID = csv_file.split('_')[-1].split('.')[0]


for i, (ra, dec) in enumerate(zip(crpix1, crpix2)):

    #this is how I will update the dictionary
    template_header['csys']['direction0']['crpix'] = numpy.array([ra, dec])

    if len(axis) == 4:
        fourth_axis = axis[3]
        if fourth_axis == 'Frequency':
            number_of_frequency = fitsheader['shape'][3]
            template_header['shap'] = numpy.array([2048, 2048, 1, number_of_frequency])

        third_axis = axis[2]
        if third_axis == 'Frequency':
            number_of_frequency = fitsheader['shape'][2]
            template_header['shap'] = numpy.array([2048, 2048, number_of_frequency, 1])

    if len(axis) == 3:
        third_axis = axis[2]
        if third_axis == 'Frequency':
            number_of_frequency = fitsheader['shape'][2]
            template_header['shap'] = numpy.array([2048, 2048, number_of_frequency])
        else:
            template_header['shap'] = numpy.array([2048, 2048, 1])

    outputname = outdir + 'tile-%s-%d.image'%(SB_ID, pixel_ID[i])

    # tiling
    imregrid(imagename=imagecasa, template=template_header, output=outputname, axes=[0, 1],  interpolation='cubic')
    #convert .image to fits files
    exportfits(imagename=outputname, fitsimage=outputname.split('.image')[0]+ '.fits')
    # delete all fits files.
    os.system('rm -rf %s'%outputname)

os.system('rm -rf %s'%imagecasa)

