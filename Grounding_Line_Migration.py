import IceSatHDF5Unpacker as unpack
import FlowUnpacker as flow
import GLineUnpacker as gline
import Visualizations as visualize
import geopandas as gpd

path = r'C:\Users\Nacnu\Documents\IceSat2\backgrounddata\ATA_adm0.shp'

#flowdatabase = flow.Database("data/antarctic_ice_vel_phase_map_v01.nc")
#flowdatabase.compute_all()

#dataset = unpack.Dataset(r"C:\Users\Nacnu\Documents\IceSat2")

#granule = dataset.openfilename("ATL06_20181114050603_07110110_005_01.h5")
#lat, lon, time, dh_fit_dx, dh_fit_dx_sigma, metadata, granuledata = granule.laser1l.getTrackData()
#dh_fit_dy = granule.laser1l.returnAcrossTrackSlope()

#xyindices, trackx, tracky = unpack.Basemap.latlonindex_to_grid(lat, lon, flowdatabase)
#flow_angle, flow_angle_error = flowdatabase.rearrange_angle_data(xyindices)
#azumith_in_xy = unpack.Basemap.angleTransform(trackx, tracky)
#visualize.Angles.visualize_ang_xy(gpd.read_file(path), trackx, tracky, azumith_in_xy, 1000, 100) 
#flowslopes = flowdatabase.get_flow_slopes(dh_fit_dx, dh_fit_dy, dh_fit_dx_sigma, azumith_in_xy, flow_angle, flow_angle_error)
#visualize.Angles.flowslope(dh_fit_dx[2000], dh_fit_dy[2000], azumith_in_xy[2000], flow_angle[2000], flowslopes[2000])

glinedata = gline.GLineData("grounding_line_data/InSAR_GL_Antarctica_v02.shp")
visualize.Polygons.line_buffer(gpd.read_file(path), glinedata.linedata, glinedata.linepolygon)