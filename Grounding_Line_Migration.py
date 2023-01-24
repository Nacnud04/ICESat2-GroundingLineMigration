import IceSatHDF5Unpacker as unpack
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import polygon as Polygon
import time

#path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'
#basemap = unpack.Basemap(path)

starttime = time.time()

dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

midtime = time.time()

granule = dataset.openfilename("ATL06_20181114050603_07110110_005_01.h5")
coverage = granule.captureCoverage()

#trackdata = dataset.returnTrackDataByPolygon(coverage)
print(list(dataset.returnFilenamesByPolygon(coverage)))
print(f"Total runtime: {time.time()-starttime}")
print(f"Midpoint time: {time.time()-midtime}")

#basemap.plotTracks(list(trackdata), datatype=tuple)