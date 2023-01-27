import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import pandas as pd
import geopandas as gpd
import numpy as np
import math
from cartopy import crs as ccrs

#path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'
#basemap = unpack.Basemap(path)

flowdatabase = flow.Database("data/antarctic_ice_vel_phase_map_v01.nc")
flowdatabase.compute_angle()

dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

granule = dataset.openfilename("ATL06_20181114050603_07110110_005_01.h5")
lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = granule.laser1l.getTrackData()

def latlonindex_to_grid(lat, lon, time, dh_fit_dh, dh_fit_dx_sigma):
    crs = ccrs.SouthPolarStereo()
    crs_proj4 = crs.proj4_init

    track_df = pd.DataFrame({'time': time, 'slope': dh_fit_dx, 'slope_error': dh_fit_dx_sigma, 'lat': lat, 'lon':lon})
    track_gdf = gpd.GeoDataFrame(track_df, geometry=gpd.points_from_xy(track_df.lon, track_df.lat), crs="EPSG:4326")

    basemap_gpd = track_gdf.to_crs(crs_proj4)
    #basemap_gpd.geometry.to_numpy()

    # find most local ice velocity vector
    def find_nearest(array,value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
            return idx
        else:
            return idx

    def find_nearest_grid(arrx, arry, valx, valy):
        idx = find_nearest(arrx, valx)
        idy = find_nearest(arry, valy) 
        return idx, idy

    xyindices = basemap_gpd.apply(lambda row: find_nearest_grid(flowdatabase.x, np.flip(flowdatabase.y), row.geometry.x, row.geometry.y), axis=1)
    return xyindices

def rearrange_flow_data(flowdatabase, xyindices):
    ERR = np.array([[flowdatabase.ERRX[idx, idy], flowdatabase.ERRY[idx, idy]] for idx, idy in xyindices])
    VEL = np.array([[flowdatabase.VX[idx, idy], flowdatabase.VY[idx, idy]] for idx, idy in xyindices])
    if flowdatabase.angle:
        ANG = np.array([flowdatabase.angle[idx, idy] for idx, idy in xyindices])
    if flowdatabase.speed:
        SPD = np.array([flowdatabase.speed[idx, idy] for idx, idy in xyindices])
    if flowdatabase.error:
        TOTERR = np.array([flowdatabase.error[idx, idy] for idx, idy in xyindices])
    if flowdatabase.angle_error:
        ANGERR = np.array([flowdatabase.angle_error[idx, idy] for idx, idy in xyindices])
