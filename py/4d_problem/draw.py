#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python3.9
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import draw_ndf

script, filename = sys.argv
draw_ndf.main(filename)
#optimized_params = 4
#magnet_dim = len(pd.read_csv(filename).columns)-optimized_params
#qNom = np.zeros(magnet_dim) 
#qNom_orig = qNom
#quads = []
#for i in range(magnet_dim):
#    quads.append("y{}".format(i))
#columns = quads
#columns.append("resol_FP2")
#columns.append("e_xangle_FP2")
#columns.append("resol_FP3")
#columns.append("e_xangle_FP3")
#
#df = pd.read_csv(filename,names=columns)
#print(df)
#df["x"] = 0
#df["avg_y0"] = 0
#df["avg_y1"] = 0
#df["resol_max"] = 1000000
#df["xwidth_e_min"] = 1000000
#df["xangle_e_min"] = 1000000
#df["xangle_xwidth_min"] = 1000000
##df["resol"] = df['resol'] * 1.0 / 4.419411469795324
##df['xangle_xwidth'] = df['xangle_xwidth'] * 1.0 / 2.7701491204695756
#max_y0, max_y1 = 0,0
#max_y = np.zeros(magnet_dim)
#n_islands = 1
#pop_size = 84 
##pop_size = 100 
#pop_gen = pop_size*n_islands
#for i in range(len(df["x"])):
#    df.loc[i,"x"] = int((i-pop_gen)//pop_gen)+1 
##    avg0 = np.average(df.loc[df["x"]==df["x"][i],"y0"])
##    avg1 = np.average(df.loc[df["x"]==df["x"][i],"y1"])
##    df.loc[df["x"]==df["x"][i],"avg_y0"]=avg0       
##    df.loc[df["x"]==df["x"][i],"avg_y1"]=avg1       
##    print(df["resol_max"].min())
#    df.loc[i,"resol_max"] = min(df["resol_max"].min(),df.loc[i,"resol_FP2"].min())
##    df.loc[i,"xwidth_e_min"] = min(df["xwidth_e_min"].min(),df.loc[i,"xwidth_e"].min())
#    df.loc[i,"xangle_e_min"] = min(df["xangle_e_min"].min(),df.loc[i,"e_xangle_FP2"].min())
##    df.loc[i,"xangle_xwidth_min"] = min(df["xangle_xwidth_min"].min(),df.loc[i,"xangle_xwidth"].min())
#    if df.loc[i,"resol_FP2"].min() == df["resol_max"].min():
#        max_y = df.iloc[i][0:magnet_dim]
#print(df)
#print(np.power(2,max_y))
#print(max_y)
#fig, axs = plt.subplots(4,4)
#plot_x, plot_y = 0,0
#axs[plot_y,plot_x].plot(df["x"],df["resol_max"],color="red",linestyle="dashed")
#axs[plot_y,plot_x].set_title('best resolution')
#axs[plot_y,plot_x].axes.set_yscale("log")
#plot_x += 1
##axs[plot_y,plot_x].plot(df["x"],df["xwidth_e_min"],color="blue",linestyle="dashed")
##axs[plot_y,plot_x].set_title('best xwidth_e')
##axs[plot_y,plot_x].axes.set_yscale("log")
##plot_x += 1
#axs[plot_y,plot_x].plot(df["x"],df["xangle_e_min"],color="green",linestyle="dashed")
#axs[plot_y,plot_x].set_title('best xangle_e')
#axs[plot_y,plot_x].axes.set_yscale("log")
##plot_x, plot_y = 0, plot_y+1 
##axs[plot_y,plot_x].plot(df["x"],df["xangle_xwidth_min"],color="purple",linestyle="dashed")
##axs[plot_y,plot_x].set_title('best xangle_xwidth')
##axs[plot_y,plot_x].axes.set_yscale("log")
#plot_x += 1
#for i in range(magnet_dim):
##    axs[0,1].plot(df["x"],df["y{0}".format(i)],marker="x",color="black",linestyle="")
##plt.plot(df["x"],df["avg_y"],color="black")
##plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773,color="blue",linestyle="dashed")
##plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773*1.5,color="red",linestyle="dashed")
##plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773*0.5,color="red",linestyle="dashed")
##plt.show()
#    if plot_x > 3:
#        plot_x, plot_y = 0, plot_y+1 
#    axs[plot_y,plot_x] = df['y{0}'.format(i)].plot.hist(ax=axs[plot_y,plot_x],bins=50)
#    axs[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
#    axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
#    axs[plot_y,plot_x].axes.yaxis.set_visible(False)
#    axs[plot_y,plot_x].set_title("q{0}".format(i+1))
#    plot_x += 1
#
##df2 = pd.read_csv('test40_0.5_0.5_gaussian_0.15.csv')
##axs[plot_y,plot_x] = df.loc[df["resol"]<1e5]['resol'].plot.hist(ax=axs[plot_y,plot_x],bins=20)
##axs[plot_y,plot_x].axes.yaxis.set_label("resol")
##axs[plot_y,plot_x].axes.set_yscale("log")
##axs[plot_y,plot_x].axes.set_xscale("log")
##axs[plot_y,plot_x].set_title("hist of resol")
##print(np.average(df['resol']))
#fig.tight_layout()
#plt.savefig(filename + ".png")



