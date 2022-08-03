# SkyTiles

#### Purpose : Extract tiles for a specific SB observation.

#### Definitions of Content: 

	PyMapSkyTiles.py : Used for extracting tile locations and IDs for a specific SB, and outputs the csv file 
	                     needed for tiling. The csv file contains CRPIX1/2 values.  

	PyTiling.py      : The actual script for generating the tile images and the corresponding header files.   
	
	Config_SkyMapTiles.json: Input to PyMapSkyTiles.py.  
	
	Config_Tiling.json : Input to PyTiling.py.   
	
	NB: With these Json files, you are able to specify the path and filenames to i) footprints, ii)   
	images, and iii) csv files, also provide tile parameters such as size (naxis), image sample size (cdelt),  
	Healpix nside, etc.   

#### The scripts are executed as follows:
	
	
							#./PyMapSkyTiles.py -j Config_SkyMapTiles.py
							
							#./PyTiling.py -j Config_Tiling.py
							

#### Config_SkyMapTiles.json input Json Definitions:

	path_footprints   : Path to footprints.
    run_all_footprints: true/false. If true, look inside path_footprints and extract Tile information for all footprints. 
    footprint_files   : If run_all_footprints is false, then specify a specific footprint file(s). 
    HPX_nside         : Healpix NSIDE 
    tile_naxis        : Tile size (tile_naxis by tile_naxis).
    tile_cdelt        : Tile pixel size.
    beam_radius       : Size of the observation's primary beam.
    beam_sample_points: Sampling the circumference of the beam. Used for determining pixel locations 
	                    (ensures that no tile ID is missed).
    number_of_beams   : The number of beams per SB. Most relevant for phased arrays. 
    outfile_prefix    : The prefix name to use for the output csvs files.
    generate_ds9regions: if true, generate a region file for the SB.
     
	Output of executing PyMapSkyTiles.py: 
		i) csv file for each SB containing tile IDs, CRPIX1/2 and CRVAL1/2.The output name takes 
		the form: 'output_prefix_TileConfig_SBID.csv'
        ii) one csv file containing all tile IDs corresponding to the input footprints. This file 
		is important to note which tiles consists of contribution from multiple SBs.  The output 
		name takes the form: 'output_prefix_TileRepeat_SBID.csv'.
									  

#### Config_Tiling.json input Json Definitions:

	input_image: Information about the input image. Insert path and name (below).
         path: 
         name: 
		 
    input_tile_config: CSV file containing tile ID, CRPIX1/2 information. Specify the path and name to these files (below). 
	            Note: The input csv file must take the format: 'filename_SBID.csv'. We use this file to also extract SBID. 
        path:
        SB_tile:
   
    output_tile: Output tile images and name. 
        path : Path to save the tiles and headers.
        naxis: Tile size should be the same as tile_naxis in Config_SkyMapTiles.json.
        cdelt: Tile pixel sample. Should be the same as tile_cdelt in Config_SkyMapTiles.json.
        fits_prefix: The prefix to give the resulting tile images and headers.         
			The output tiles/headers will take the form: 
			'fits_prefix-stokes-SBID-header-TileID.hdr' and 'fits_prefix-stokes-SBID-tile-TileID.fits.     
			The stokes parameter is obtained from the input fits header with crval3/4: 1 is I, 
			2 is Q, 3 is U and 4 is V.
					 
					 

