import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

script, filename = sys.argv

qNom = np.array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )

#df = pd.read_csv(filename,names=["y0","y1","resol"])
df = pd.read_csv(filename,names=["y0","y1","y2","y3","y4","y5","y6","resol"])
print(df)
df["x"] = 0
df["avg_y0"] = 0
df["avg_y1"] = 0
df["resol_max"] = 0
max_y0, max_y1 = 0,0
n_islands = 1
pop_size = 100
pop_gen = pop_size*n_islands
for i in range(len(df["x"])):
    df.loc[i,"x"] = int((i-pop_gen)//pop_gen)+1 
#    avg0 = np.average(df.loc[df["x"]==df["x"][i],"y0"])
#    avg1 = np.average(df.loc[df["x"]==df["x"][i],"y1"])
#    df.loc[df["x"]==df["x"][i],"avg_y0"]=avg0       
#    df.loc[df["x"]==df["x"][i],"avg_y1"]=avg1       
#    print(df["resol_max"].min())
    df.loc[i,"resol_max"] = min(df["resol_max"].min(),df.loc[i,"resol"].min())
    if df.loc[i,"resol"].min() == df["resol_max"].min():
        max_y = df.iloc[i][0:7]
print(df)
print(max_y)
fig, axs = plt.subplots(3,3)
plot_x, plot_y = 0,0
axs[plot_x,plot_y].plot(df["x"],df["resol_max"],color="red",linestyle="dashed")
plot_num = 1
plot_x += 1
for i in range(7):
#    axs[0,1].plot(df["x"],df["y{0}".format(i)],marker="x",color="black",linestyle="")
#plt.plot(df["x"],df["avg_y"],color="black")
#plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773,color="blue",linestyle="dashed")
#plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773*1.5,color="red",linestyle="dashed")
#plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773*0.5,color="red",linestyle="dashed")
#plt.show()
    if plot_x > 2:
        plot_x, plot_y = 0, plot_y+1 
    axs[plot_y,plot_x] = df['y{0}'.format(i)].plot.hist(ax=axs[plot_y,plot_x],bins=50)
    axs[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
    axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
    axs[plot_y,plot_x].axes.yaxis.set_visible(False)
    axs[plot_y,plot_x].set_title("q{0}".format(i+1))
    plot_x += 1

#df2 = pd.read_csv('test40_0.5_0.5_gaussian_0.15.csv')
axs[plot_y,plot_x] = df['resol'].plot.hist(ax=axs[plot_y,plot_x],bins=20)
axs[plot_y,plot_x].axes.yaxis.set_label("resol")
axs[plot_y,plot_x].axes.set_yscale("log")
print(np.average(df['resol']))
fig.tight_layout()
plt.savefig(filename + ".png")



