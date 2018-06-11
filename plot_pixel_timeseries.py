import datetime
import numpy as np
import matplotlib.pyplot as plt

import gdal

from pathlib import Path


def do_plot(the_path, ax, params = None, marker = '-'):
    path=Path(the_path)
    #params=["w_nir", "x_nir", "a_nir",
    #    "w_vis", "x_vis", "a_vis"]
    if params == None:
        params=["TeLAI", "w_vis", "a_vis"]
    
    labels = {"TeLAI":"Transformed LAI", "w_nir":"leaf single scattering albedo in NIR", 
              "x_nir":"leaf assymetry factor in NIR", "a_nir":"background albedo in NIR", 
              "w_vis":"leaf single scattering albedo in NIR", "x_vis":"leaf assymetry factor in NIR", 
              "a_vis":"background albedo in vis"}

    data = {}
    uncs = {}
    dates =[]
    for param in params:
        ll = []
        uu = []
        
        for unc in [True, False]:
            if unc:
                files = np.sort([f for f in path.glob(f"{param}_A2017???_unc.tif")])
            else:
                files = np.sort([f for f in path.glob(f"{param}_A2017???.tif")])
            for fich in files:
                g = gdal.Open(fich.as_posix())
                if not unc and param == "TeLAI":
                    dates.append(int(datetime.datetime.strptime(fich.name.split("_")[1][1:],"%Y%j.tif").strftime("%j")))
                elif not unc and param == params[0]:
                    dates.append(int(datetime.datetime.strptime(fich.name.split("_")[2][1:],"%Y%j.tif").strftime("%j")))
                d = g.ReadAsArray()[698, 1230]
                if unc:
                    uu.append(d)
                else:
                    ll.append(d)
        ax.plot(np.array(dates), np.array(ll), marker, label=labels[param])
        ax.fill_between(np.array(dates), np.array(ll) - np.array(uu), np.array(ll)+np.array(uu),
                            color="0.8", alpha=0.5)
        ax.set_ylim([-0.1, 1.1])
        data[param] = np.array(ll)
        uncs[param] = np.array(uu)
        
        
    ax.legend(loc="best")
    return data, uncs, dates
    
'''    
fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, sharey=True, figsize=(8,20))
axs = axs.flatten()
paths = ["/home/npounder/output/kafka/28Mar/kafkaout_no_prop_HessCorr",
        "/home/npounder/output/kafka/28Mar/kafkaout_prop_LAI_Q0-1/",
        "/home/npounder/output/kafka/28Mar/kafkaout_prop_LAI_Q0-25/",
        "/home/npounder/output/kafka/28Mar/kafkaout_prop_LAI_Q0-5/",
        "/home/npounder/output/kafka/28Mar/kafkaout_prop_LAI_Q0-75"]
titles = ["no propagation", "0.1", "0.25", "0.5", "75"]
for ax, the_path, title in zip(axs, paths, titles):
    do_plot(the_path, ax)
    ax.set_title(title)

'''

