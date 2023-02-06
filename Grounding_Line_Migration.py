import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import GLineUnpacker as gline
import Visualizations as visualize
import geopandas as gpd
import numpy as np

dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

groundline = gline.GLineData("grounding_line_data/InSAR_GL_Antarctica_v02.shp")
polygons = groundline.linepolygon.set_crs(crs="+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

files = np.array(dataset.returnFilenamesByPolygon(polygons))
files = np.unique(files)

print(f"Original filecount: {len(dataset.datafiles)}  |  New filecount: {len(files)}")
