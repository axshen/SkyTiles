import pyCASATILE


def test_casa_tile_sb_10040():
    """Test running on AusSRC Carnaby with expected file system.

    """
    pyCASATILE.main([
        '-i', '/mnt/shared/home/ashen/POSSUM/data/3chan/image.restored.i.SB10040.contcube_3chan.fits',
        '-m', '/mnt/shared/home/ashen/POSSUM/development/SkyTiles/PoSSUM_2156-54.csv',
        '-o', '/mnt/shared/home/ashen/POSSUM/development/SkyTiles/tiles/',
        '-j', '/mnt/shared/home/ashen/POSSUM/development/SkyTiles/Config_CASATILE.json'
    ])
