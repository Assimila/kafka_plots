import datetime as dt
import numpy as np
import gdal
import matplotlib.pyplot as plt
import pickle

from pathlib import Path
cmap_s = plt.cm.viridis
cmap_d = plt.cm.RdBu_r

def filename(filepath,parameter,date, extension='tif', unc = False):
    doy = date.timetuple().tm_yday
    year = date.year
    if unc:
        fname = "{}/{}_A{}{:0>3}_unc.{}".format(filepath,parameter,
                                                year, doy, extension)
    else:
        fname = "{}/{}_A{}{:0>3}.{}".format(filepath,parameter,year, doy, extension)
    return fname


def transform_LAI(t_LAI):
    return -2*np.log(t_LAI)


def plot_image(filepath, ax, date, vmin = -2, vmax = 6,
               xlim = None, ylim=None,parameter =" TeLAI",
               convertLAI="True", unc=False):

    fname = filename(filepath,parameter, date, extension='tif', unc=unc)
    print(fname)
    tif = gdal.Open(fname)
    if tif is None: ## need to raise an error here
        raise ValueError("Image not found {}".format(fname))
    var = tif.ReadAsArray()
    if convertLAI is True:# and parameter == "TeLAI":
        var = transform_LAI(var)
    im = ax.imshow(var, vmin=vmin, vmax=vmax, cmap=cmap_s, interpolation = 'none')
    if convertLAI and (parameter == "TeLAI" or parameter == "lai"):
        ax.set_title("LAI {}".format(date.date()))
    else:
        ax.set_title(parameter.format(date.date()))
    plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    return ax


def plot_image_tseries(filepath, axs, time_grid, vmin = -2, vmax = 6,
                           xlim = None, ylim=None, parameter ="TeLAI",
                       convertLAI=False, unc=False):
    for ax, date in zip(axs, time_grid):
        try:
            plot_image(filepath, ax, date, vmin=vmin, vmax=vmax,
                       xlim=xlim, ylim=ylim, parameter=parameter,
                       convertLAI=convertLAI, unc=unc)
        except ValueError as e:
            print(e)
            continue
    return axs


def plot_image_tseries_diff(filepath1, filepath2, axs, time_grid, vmin = -3, vmax = 3,
                           xlim = None, ylim=None):
    # This has not been updated in a while...
    for ax, date in zip(axs, time_grid):
        doy = date.timetuple().tm_yday
        fname = filename(filepath1 ,parameter, date, extension='tif', unc=unc)
        tif1 = gdal.Open(fname)
        if tif1 is None: continue
        LAI1 = -2*np.log(tif1.ReadAsArray())
        fname = filename(filepath2 ,parameter, date, extension='tif', unc=unc)
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


def plot_scatter_tseries(filepath1, filepath2, axs, time_grid,parameter =" TeLAI",
                         convertLAI="True", unc=False, vmin=-1, vmax=6,
                         xlim = None, ylim=None):
    for ax, date in zip(axs, time_grid):
        fname = filename(filepath1,parameter, date, extension='tif', unc=unc)
        print(fname)
        tif1 = gdal.Open(fname)
        if tif1 is None: continue
        fname = filename(filepath2,parameter, date, extension='tif', unc=unc)
        print(fname)
        tif2 = gdal.Open(fname)
        if tif2 is None: continue
        a1 = tif1.ReadAsArray()
        a2 = tif2.ReadAsArray()
        if xlim is not None and ylim is not None:
            a1 = a1[ylim[0]:ylim[1], xlim[0]:xlim[1]]
            a2 = a2[ylim[0]:ylim[1], xlim[0]:xlim[1]]
        a1 = a1.flatten()
        a2 = a2.flatten()
        if convertLAI is True:
            a1 = transform_LAI(a1)
            a2 = transform_LAI(a2)
        passer1 = (a1<100) * (a1>-100)
        passer2 = (a2<100) * (a2>-100)
        #passer2 = np.isfinite(y_data[doy])
        passer = passer1*passer2
        im = ax.hexbin(-2*np.log(a1[passer]),-2*np.log(a2[passer]),
                       gridsize=50, bins="log", mincnt=5,cmap=cmap_s)
        ax.plot([vmin, vmax], [vmin, vmax])
        ax.set_title(date)
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)
        plt.colorbar(im, ax=ax,fraction=0.046, pad=0.04)
        tif1 = None
        tif2 = None
    return axs

def plot_image_diff(filepath1, filepath2, date, ax, vmin = -2, vmax = 6,
                           xlim = None, ylim=None, parameter ="LAI"):
    # This has not been updated in a while...
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


def plot_pixel_tseries(ax, data, unc, dates, parameter ="TeLAI",
                       marker = '-', convertLAI = False, plot_unc=True, ):
    labels = {"TeLAI":"Transformed LAI",
              "lai":"Transformed LAI",
              "w_nir":"leaf single scattering albedo in NIR",
              "x_nir":"leaf asymmetry factor in NIR",
              "a_nir":"background albedo in NIR",
              "w_vis":"leaf single scattering albedo in NIR",
              "x_vis":"leaf asymmetry factor in NIR",
              "a_vis":"background albedo in vis"}
    l_unc = data-unc
    u_unc = data+unc
    if convertLAI is True:
        data = -2*np.log(data)
        l_unc = -2.*np.log(l_unc)
        u_unc = -2.*np.log(u_unc)
        l_unc[np.isnan(l_unc)] = 10
        labels["TeLAI"] = "LAI"
        labels["lai"] = "LAI"
    p = ax.plot(np.array(dates), np.array(data), marker, label=labels[parameter])

    if plot_unc:
        #color = "0.8"
        color = (p[0].get_color())
        ax.set_title(labels[parameter])
        ax.fill_between(np.array(dates), l_unc, u_unc,
                    color=color, alpha=0.4)
    return ax


def plot_pixel_tseries_s2(ax, data, unc, dates, parameter="lai",
                          marker='-', convertLAI=False, error="shading",
                          line='None', label='Auto'):
    labels = {"lai": "Transformed LAI",
              'n': "n",
              'cab': "cab",
              'car': "car",
              'cbrown': "cbrown",
              'cw': "cw",
              'cm': "cm",
              'ala': "ala",
              'bsoil': 'bsoil',
              'psoil': "psoil"}
    # try:
    # for d,u, l,t, td  in zip(data, unc, (data-unc), -2.*np.log(data-unc), -2.*np.log(data)):
    #    if not np.isnan(t):
    #        print (d,u,l,t, td)
    l_unc = data + np.sqrt(unc)
    # except ValueError:
    #    print ("waring - unc calc error {}".format(unc))
    #    l_unc = data
    # try:
    u_unc = data - np.sqrt(unc)
    # except ValueError:
    #    print ("waring - unc calc error {}".format(unc))
    #    u_unc = data
    if convertLAI is True:
        data = -2 * np.log(data)
        l_unc = -2. * np.log(l_unc)
        u_unc = -2. * np.log(u_unc)
        l_unc[np.isnan(l_unc)] = 10
        labels["lai"] = "LAI"
    maskprior = np.where(abs(data - 4.0) > 0.0000001)
    if label == 'Auto': label = labels[parameter]
    if error == "bar":
        uncs = np.array([data[maskprior] - l_unc[maskprior], u_unc[maskprior] - data[maskprior]])
        p = ax.errorbar(np.array(np.array(dates)[maskprior]), np.array(data)[maskprior],
                        yerr=uncs, linestyle=line,
                        marker=marker, label=label)  # , color = 'k')
    else:
        p = ax.plot(np.array(np.array(dates)[maskprior]), np.array(data)[maskprior],
                    marker=marker, linestyle=line, label=label)
    if error == "shading":
        color = (p[0].get_color())
        ax.fill_between(np.array(dates), l_unc, u_unc,
                        color=color, alpha=0.4, label='_nolegend_')

    ax.set_title(labels[parameter])
    return ax


def extract_pixel(filepath, year, x,y, params=None, outfile=None):
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
    doys = []
    for param in params:
        tseries = []
        uncertainty = []
        
        files = np.sort([f for f in path.glob(f"{param}_A{year}???.tif")])
        for f in files:
            g = gdal.Open(f.as_posix())
            if param == params[0]:
                # Below gives day of year as an int.
                doys.append(int(dt.datetime.strptime(f.name.split("_")[len(param.split("_"))][1:],"%Y%j.tif").strftime("%j")))
                # Below gives a datetime object
                dates.append(dt.datetime.strptime(f.name.split("_")[len(param.split("_"))][1:],"%Y%j.tif"))
                #print(dates)
                #print(len(param.split("_")))
            d = g.ReadAsArray()[y,x]
            tseries.append(d)

        files = np.sort([f for f in path.glob(f"{param}_A{year}???_unc.tif")])
        for f in files:
            g = gdal.Open(f.as_posix())
            d = g.ReadAsArray()[y,x]
            uncertainty.append(d)
        
        data[param] = np.array(tseries)
        uncs[param] = np.array(uncertainty)
        
        pickle.dump((data, uncs, dates, doys), open(outfile, "wb"))

    return data, uncs, dates, doys
    

def get_pixel(filepath, year, x, y, params = None, recreate_file = False):
    """
    Get retrieved parameters and uncertainties for a single pixel
    This will open the data from a pickle file if it exists, otherwise
    it calls extract_pixel, which will get the data and also store a pickle file
    for next time.
    """
    if params == None:
        params = ['n', 'cab', 'car', 'cbrown', 'cw', 'cm',
                  'lai', 'ala', 'bsoil', 'psoil']
    file = filepath+"/pixel_{}_{}.pkl".format(x,y)
    try:
        data, uncs, dates, doys = pickle.load(open(file, 'rb'))
    except (FileNotFoundError, ValueError):
        recreate_file = True
    if recreate_file:
        print('creating '+file)
        data, uncs, dates, doys = extract_pixel(filepath, year, x, y, params=params, outfile=file)
    return data, uncs, dates, doys

