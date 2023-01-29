import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import numpy as np
import math

#path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'
#basemap = unpack.Basemap(path)

flowdatabase = flow.Database("data/antarctic_ice_vel_phase_map_v01.nc")
flowdatabase.compute_all()

dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

granule = dataset.openfilename("ATL06_20181114050603_07110110_005_01.h5")
lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = granule.laser1l.getTrackData()
azumith = granule.laser1l.returnAzimuth()
dh_fit_dy = granule.laser1l.returnAcrossTrackSlope()

xyindices = unpack.Basemap.latlonindex_to_grid(lat, lon, flowdatabase)
flow_angle, flow_angle_error = flowdatabase.rearrange_angle_data(xyindices)
azumith_in_xy = unpack.Basemap.angleTransform(lon, lat, azumith)
flowslopes = flowdatabase.get_flow_slopes(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, azumith_in_xy, flow_angle, flow_angle_error)
print(flowslopes)