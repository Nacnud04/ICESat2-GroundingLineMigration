# Visualizations file
# Used to visualize data
from math import sin, cos, pi
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import geopandas as gpd
import numpy as np

class Angles:

    def __init__(self):
        pass

    @staticmethod
    def visualize_ang_xy(basemap, xs, ys, azumiths, vectorlength:float, fill:int):
        crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

        basemap_gpd = basemap.to_crs(crs_proj4)

        maxx, minx = max(xs)*1.25, min(xs)/1.25
        maxy, miny = max(ys)*1.25, min(ys)/1.25

        ax = basemap_gpd.plot(color="white", edgecolor="black")
        xs = [x for i, x in enumerate(xs) if i % fill == 0]
        ys = [y for i, y in enumerate(ys) if i % fill == 0]
        azumiths = [a for i, a in enumerate(azumiths) if i % fill == 0]
        ax.scatter(xs, ys)
        for i in range(len(xs)):
            ax.plot([xs[i], xs[i]+vectorlength*cos(azumiths[i])], [ys[i], ys[i]+vectorlength*sin(azumiths[i])], color = "red")
        plt.xlim(minx, maxx)
        plt.ylim(miny, maxy)
        plt.show()

    @staticmethod 
    def flowslope(dh_fit_dx, dh_fit_dy, azumith, flowdirec, flowslope):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        trackslope = np.array([cos(azumith), sin(azumith), dh_fit_dx])
        ts_hat = trackslope / np.linalg.norm(trackslope)
        acrosstrackslope = np.array([sin(azumith), -1*cos(azumith), dh_fit_dy])
        ats_hat = acrosstrackslope / np.linalg.norm(acrosstrackslope)
        flowslope = np.array([cos(flowdirec), sin(flowdirec), flowslope])
        f_hat = flowslope / np.linalg.norm(flowslope)
        ax.plot([0, cos(azumith)], [0, sin(azumith)], [0, 0], color="red")
        ax.plot([0, ts_hat[0]], [0, ts_hat[1]], [0, ts_hat[2]], color="red")
        ax.plot([0, ats_hat[0]], [0, ats_hat[1]], [0, ats_hat[2]], color="green")
        ax.plot([0, cos(flowdirec)], [0, sin(flowdirec)], [0, 0], color="blue")
        ax.plot([0, f_hat[0]], [0, f_hat[1]], [0, f_hat[2]], color="blue")
        plt.show()


class Polygons():

    def __init__(self):
        pass

    @staticmethod 
    def line_buffer(basemap, line, buff):
        crs_proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
        basemap_gpd = basemap.to_crs(crs_proj4)
        fig, ax = plt.subplots(1,1,sharex=True,sharey=True,figsize=(11,11))
        basemap_gpd.plot(ax = ax, color="white", edgecolor="black")
        buff.plot(ax = ax, color = "lightblue", edgecolor = "steelblue", alpha = 0.5)
        line.plot(ax = ax, color = "crimson")
        plt.show()