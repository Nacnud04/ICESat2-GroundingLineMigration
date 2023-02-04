import geopandas as gpd

class GLineData():

    def __init__(self, filepath):
        
        data = gpd.read_file(filepath)
        self.crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
        self.linedata = data.to_crs(self.crs_proj4)
        
        self.linepolygon = self.generateBuffer(100000)

    def generateBuffer(self, buffersize):

        linepolygon = self.linedata.buffer(buffersize)
        linepolygon = gpd.GeoDataFrame(
              geometry=[linepolygon.unary_union]).explode(
              index_parts=False).reset_index(
              drop=True)
        return linepolygon