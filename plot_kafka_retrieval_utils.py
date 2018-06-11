import datetime as dt
import numpy as np
import gdal
import matplotlib.pyplot as plt
import pickle

from pathlib import Path
cmap_s = plt.cm.viridis
cmap_d = plt.cm.RdBu_r

def plot_image_tseries(filepath, axs, time_grid, vmin = -2, vmax = 6,
                           xlim = None, ylim=None, parameter ="TeLAI", convertLAI=False, unc=False):
    for ax, date in zip(axs, time_grid):
        doy = date.timetuple().tm_yday
        year = date.year
        if unc:
            fname = "{}/{}_A{}{:0>3}_unc.tif".format(filepath,parameter,year, doy)
            convertLAI = False # There is no code to correctly convert LAI uncertainty yet
        else:
            fname = "{}/{}_A{}{:0>3}.tif".format(filepath,parameter,year, doy)
        print(fname)
        tif = gdal.Open(fname)
        if tif is None: 
            print("No image")
            continue
        var = tif.ReadAsArray()
        if convertLAI:# is True and parameter == "TeLAI":
            var = -2*np.log(var)
            #var = np.exp(-0.5*var)
        im = ax.imshow(var, vmin=vmin, vmax=vmax, cmap=cmap_s, interpolation = 'none')
        ax.set_title(date.date())
        plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        tif = None
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
    return axs

def plot_image_and_unc_tseries(filepath, axs, time_grid, vmin = -2, vmax = 6,
                           xlim = None, ylim=None, parameter ="TeLAI", convertLAI=True, unc=False):
    for ax, date in zip(axs, time_grid):
        doy = date.timetuple().tm_yday
        year = date.year
        if unc:
            fname = "{}/{}_A{}{:0>3}_unc.tif".format(filepath,parameter,year, doy)
            convertLAI = False # There is no code to correctly convert LAI uncertainty yet
        else:
            fname = "{}/{}_A{}{:0>3}.tif".format(filepath,parameter,year, doy)
        print(fname)
        tif = gdal.Open(fname)
        if tif is None: 
            print("No image")
            continue
        var = tif.ReadAsArray()
        if convertLAI is True:# and parameter == "TeLAI":
            var = -2*np.log(var)
            #var = np.exp(-0.5*var)https://www.longtallsally.com/point-zero/c
        im = ax.imshow(var, vmin=vmin, vmax=vmax, cmap=cmap_s, interpolation = 'none')
        ax.set_title(date)
        plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        tif = None
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
    return axs


def plot_image_tseries_diff(filepath1, filepath2, axs, time_grid, vmin = -3, vmax = 3,
                           xlim = None, ylim=None):
    for ax, date in zip(axs, time_grid):
        doy = date.timetuple().tm_yday
        fname = "{}/TeLAI_A2016{:0>3}.tif".format(filepath1, doy)
        #print fname
        tif1 = gdal.Open(fname)
        if tif1 is None: continue
        LAI1 = -2*np.log(tif1.ReadAsArray())
        fname = "{}/TeLAI_A2016{:0>3}.tif".format(filepath2, doy)
        #print fname
        tif2 = gdal.Open(fname)
        if tif2 is None: continue
        LAI2 = -2*np.log(tif2.ReadAsArray())
        im = ax.imshow(LAI1-LAI2, vmin=vmin, vmax=vmax, cmap=cmap_d, interpolation = 'none')
        ax.set_title(date.date())
        plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        tif1 = None
        tif2 = None
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
    return axs


def plot_scatter_tseries(filepath1, filepath2, axs, time_grid, vmin = -1, vmax = 1.2):
    for ax, date in zip(axs, time_grid):
        doy = date.timetuple().tm_yday
        fname = "{}/TeLAI_A2016{:0>3}.tif".format(filepath1, doy)
        #print fname
        tif1 = gdal.Open(fname)
        if tif1 == None: continue
        fname = "{}/TeLAI_A2016{:0>3}.tif".format(filepath2, doy)
        #print fname
        tif2 = gdal.Open(fname)
        if tif2 == None: continue
        a1 = tif1.ReadAsArray().flatten()
        a2 = tif2.ReadAsArray().flatten()
        passer1 = (a1<100) * (a1>-100)
        passer2 = (a2<100) * (a2>-100)
        #passer2 = np.isfinite(y_data[doy])
        passer = passer1*passer2
        im = ax.hexbin(-2*np.log(a1[passer]),-2*np.log(a2[passer]),
                       gridsize=50, bins="log", mincnt=5,cmap=cmap_s)
        ax.set_title(date)
        ax.set_xlim(-6, 6)
        ax.set_ylim(-6, 6)
        plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        tif1 = None
        tif2 = None
    return axs

def plot_image_diff(filepath1, filepath2, date, ax, vmin = -2, vmax = 6,
                           xlim = None, ylim=None, parameter ="LAI"):
    convertLAI = False
    if parameter == "LAI":
        convertLAI = True
        parameter = "TeLAI"
    doy = date.timetuple().tm_yday
    fname = "{}/{}_A2016{:0>3}.tif".format(filepath1,parameter, doy)
    #print fname
    tif1 = gdal.Open(fname)
    var1 = tif1.ReadAsArray()
    fname = "{}/{}_A2016{:0>3}.tif".format(filepath2,parameter, doy)
    #print fname
    tif2 = gdal.Open(fname)
    var2 = tif2.ReadAsArray()
    if convertLAI is True:
        var1 = -2*np.log(var1)
        var2 = -2*np.log(var2)
    im = ax.imshow(var2-var1, vmin=vmin, vmax=vmax, cmap=cmap_d, interpolation = 'none')
    if convertLAI:
        ax.set_title("LAI")
    else: 
        ax.set_title(parameter)
    plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    tif = None
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax

def plot_image(filepath, ax, date, xlim = None, ylim=None,
               vmin = -2, vmax = 6, parameter ="LAI", convertLAI="True"):
    doy = date.timetuple().tm_yday
    year = date.year
    fname = "{}/{}_A{}{:0>3}.tif".format(filepath,parameter,year, doy)
    #print fname
    tif = gdal.Open(fname)
    var = tif.ReadAsArray()
    if convertLAI and parameter == "TeLAI":
        print("converting")
        var = -2*np.log(var)
    im = ax.imshow(var, vmin=vmin, vmax=vmax, cmap=cmap_s, interpolation = 'none')
    if convertLAI and parameter == "TeLAI":
        ax.set_title("LAI {}".format(date.date()))
    else: 
        ax.set_title(parameter)
    plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    tif = None
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


'''
def plot_pixel(filepath, date, ax, x = None, y=None, parameter ="LAI"):
    convertLAI = False
    if parameter == "LAI":
        convertLAI = True
        parameter = "TeLAI"
    doy = date.timetuple().tm_yday
    fname = "{}/{}_A2016{:0>3}.tif".format(filepath,parameter, doy)
    #print fname
    tif = gdal.Open(fname)
    var = tif.ReadAsArray()
    if convertLAI is True:
        var = -2*np.log(var)
    
    ax.plot(var[x,y])
    
    #im = ax.imshow(var, vmin=vmin, vmax=vmax, cmap=cmap_s, interpolation = 'none')
    if convertLAI:
        ax.set_title("LAI")
    else: 
        ax.set_title(parameter)
    #plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    tif = None
'''
def plot_pixel_tseries(ax, data, unc, dates, parameter ="TeLAI",
                       marker = '-', convertLAI = False):
    labels = {"TeLAI":"Transformed LAI",
              "w_nir":"leaf single scattering albedo in NIR",
              "x_nir":"leaf assymetry factor in NIR",
              "a_nir":"background albedo in NIR",
              "w_vis":"leaf single scattering albedo in NIR",
              "x_vis":"leaf assymetry factor in NIR",
              "a_vis":"background albedo in vis"}
    l_unc = data-unc
    u_unc = data+unc
    if convertLAI is True:
        data = -2*np.log(data)
        l_unc = -2.*np.log(l_unc)
        u_unc = -2.*np.log(u_unc)
        labels["TeLAI"] = "LAI"
    ax.plot(np.array(dates), np.array(data), marker, label=labels[parameter])
    ax.fill_between(np.array(dates), l_unc, u_unc,
                    color="0.8", alpha=0.5)
    return ax





def extract_pixel(filepath, x,y, params=None, outfile=None):
    path=Path(filepath)
    if params == None:
        params=["TeLAI", "w_nir", "x_nir", "a_nir",
                "w_vis", "x_vis", "a_vis"]

    if outfile == None:
        outfile = filepath+"/pixel_{}_{}.pkl".format(x,y)
    print("saving output to {}".format(outfile))

    data = {}
    uncs = {}
    dates =[]
    for param in params:
        tseries = []
        uncertainty = []
        
        files = np.sort([f for f in path.glob(f"{param}_A2017???.tif")])
        for f in files:
            g = gdal.Open(f.as_posix())
            if param == params[0]:
                print(f.name)
                if param == "TeLAI":
                    dates.append(int(dt.datetime.strptime(f.name.split("_")[1][1:],"%Y%j.tif").strftime("%j")))
                else:
                    dates.append(int(dt.datetime.strptime(f.name.split("_")[2][1:],"%Y%j.tif").strftime("%j")))
            d = g.ReadAsArray()[x,y]
            tseries.append(d)

        files = np.sort([f for f in path.glob(f"{param}_A2017???_unc.tif")])
        for f in files:
            g = gdal.Open(f.as_posix())
            d = g.ReadAsArray()[x,y]
            uncertainty.append(d)
        
        data[param] = np.array(tseries)
        uncs[param] = np.array(uncertainty)
        
        pickle.dump((data, uncs, dates), open(outfile, "wb"))

    return data, uncs, dates
    



