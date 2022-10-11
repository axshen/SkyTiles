#!/usr/bin/env python3

"""
Code to identify which tiles to mosaic for each observation.

"""

import os
import sys
import argparse
import logging
import glob
import numpy as np
import pandas as pd


logging.basicConfig(level=logging.INFO)


def get_sbids_from_filename(filename):
    """Function to extract the unique identifier from the HPX tile fits files.

    """
    return filename.split('/')[-1].split('.')[0]


def main(argv):
    # Arguments
    parser = argparse.ArgumentParser('Determine whether new tile can be completed for a given observation.')
    parser.add_argument('-i', dest='sbid', help='Input SBID', required=True)
    parser.add_argument('-m', dest='tile_map', help='Map from HPX pixels to the required SBIDs', required=True)
    parser.add_argument('-f', dest='files', help='Directory where tiling files are stored', required=True)
    parser.add_argument('-o', dest='output', help='Output directory where tiling outputs are stored', required=True)
    args = parser.parse_args(argv)

    # Read
    tiles_df = pd.read_csv(args.tile_map)
    tiles_df = tiles_df.replace({np.nan: None})

    logging.info(f'New observation with SBID {args.sbid}')

    # Get files and previously observed sbids
    files = glob.glob(os.path.join(args.files, '*.fits'))
    previous_sbids = list(set(
        [get_sbids_from_filename(f) for f in files]
    ))
    completed_sbids = set([args.sbid] + previous_sbids)
    logging.info(f'Found HPX tiles for the previous observations: {previous_sbids}')

    # Determine which pixels can be completed
    for idx, row in tiles_df.iterrows():
        hpx_pixel_id = row.values[0]
        ids = set([v for v in row.values[1:] if v is not None])
        if ids.issubset(completed_sbids) and not ids.issubset(previous_sbids):
            logging.info(f'All observations required for HPX pixel {hpx_pixel_id} are completed.')


if __name__ == '__main__':
    main(sys.argv[1:])
