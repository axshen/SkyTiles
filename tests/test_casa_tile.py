import pyCASATILE


def test_casa_tile_sb_10040():
    pyCASATILE.main([
        '-j', '/mnt/shared/home/ashen/POSSUM/development/SkyTiles/Config_CASATILE.json'
    ])
