import IceSatHDF5Unpacker as unpack
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'
basemap = unpack.Basemap(path)

dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")
datatracks = dataset.returnTrackDataByRGT(323, num=None)
basemap.plotTracks(datatracks, datatype=tuple, local=False)