import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import pandas as pd
import geopandas as gpd
import numpy as np
import math
from cartopy import crs as ccrs
from polar_convert.constants import SOUTH
from polar_convert import polar_lonlat_to_xy

#path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'
#basemap = unpack.Basemap(path)

flowdatabase = flow.Database("data/antarctic_ice_vel_phase_map_v01.nc")
flowdatabase.compute_all()

dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

granule = dataset.openfilename("ATL06_20181114050603_07110110_005_01.h5")
lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = granule.laser1l.getTrackData()
azumith = granule.laser1l.returnAzimuth()
dh_fit_dy = granule.laser1l.returnAcrossTrackSlope()

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
    return ERR, VEL, ANG, SPD, TOTERR, ANGERR

def rearrange_angle_data(flowdatabase, xyindices):
    try:
        ANG = np.array([flowdatabase.angle[idx, idy] for idx, idy in xyindices])
        ANGERR = np.array([flowdatabase.angle_error[idx, idy] for idx, idy in xyindices])
        return ANG, ANGERR
    except:
        raise Exception("Angle data for flowdatabase does not exist")

# the magic algorithm
def angNorth_to_angXY(coord, angle):
    lon, lat = coord
    wgsvec = (lon + math.cos(angle), lat + math.sin(angle))
    true_scale_lat = 71
    re = 6378.137
    e = 0.01671
    hemisphere = SOUTH
    origxy = polar_lonlat_to_xy(lon, lat, true_scale_lat, re, e, hemisphere)
    newxy = polar_lonlat_to_xy(wgsvec[0], wgsvec[1], true_scale_lat, re, e, hemisphere)
    vecxy = [newxy[i] - origxy[i] for i in range(2)]
    angle = math.atan(vecxy[1]/vecxy[0])
    return angle

def angleTransform(lon, lat, azumith):
    angles= np.array([angNorth_to_angXY((lon[i], lat[i]), azumith[i]) for i in range(len(lon))])
    return angles

# known to work perfectly.
def calc_along_flow_slope(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, slope_azumith, flow_ang, flow_ang_err):
    slopevec = np.array([math.cos(slope_azumith), math.sin(slope_azumith), dh_fit_dx])
    acrossslopevec = np.array([math.sin(slope_azumith), -1*math.cos(slope_azumith), dh_fit_dy])
    flowvec = np.array([math.cos(flow_ang), math.sin(flow_ang), 0])
    slopeflowvec = (flowvec.dot(slopevec)/slopevec.dot(slopevec))*slopevec+(flowvec.dot(acrossslopevec)/acrossslopevec.dot(acrossslopevec))*acrossslopevec
    slopeflowgrade = slopeflowvec[2] / math.sqrt(slopeflowvec[0]**2 + slopeflowvec[1]**2)
    if dh_fit_dy >= 1e30:
        return None
    else:
        return slopeflowgrade

def get_flow_slopes(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, slope_azumith, flow_ang, flow_ang_err):
    along_flow_slope = np.array([calc_along_flow_slope(dh_fit_dx[i], dh_fit_dy[i], dh_fit_dx_sigma[i], slope_azumith[i], flow_ang[i], flow_ang_err[i]) for i in range(len(flow_ang))])
    return along_flow_slope

xyindices = unpack.Basemap.latlonindex_to_grid(lat, lon, flowdatabase)
flow_angle, flow_angle_error = rearrange_angle_data(flowdatabase, xyindices)
azumith_in_xy = angleTransform(lon, lat, azumith)
flowslopes = get_flow_slopes(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, azumith_in_xy, flow_angle, flow_angle_error)