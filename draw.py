import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

script, filename = sys.argv

qNom = np.array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )

df = pd.read_csv(filename,names=["y0","y1","resol"])
print(df)
df["x"] = 0
df["avg_y0"] = 0
df["avg_y1"] = 0
df["resol_max"] = 0
for i in range(len(df["x"])):
    df.loc[i,"x"] = int((i-8)//3)+1 
    avg0 = np.average(df.loc[df["x"]==df["x"][i],"y0"])
    avg1 = np.average(df.loc[df["x"]==df["x"][i],"y1"])
    df.loc[df["x"]==df["x"][i],"avg_y0"]=avg0       
    df.loc[df["x"]==df["x"][i],"avg_y1"]=avg1       
    print(df["resol_max"].min())
    df.loc[i,"resol_max"] = min(df["resol_max"].min(),df.loc[i,"resol"].min())

print(df)
fig, axs = plt.subplots(2,2)

axs[0,0].plot(df["x"],df["y0"],marker="x",color="black",linestyle="")
#plt.plot(df["x"],df["avg_y"],color="black")
#plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773,color="blue",linestyle="dashed")
#plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773*1.5,color="red",linestyle="dashed")
#plt.plot(df["x"],np.zeros(len(df["x"]))-0.39773*0.5,color="red",linestyle="dashed")
#plt.show()

axs[1,0].plot(df["x"],df["resol_max"],color="red",linestyle="dashed")
axs[0,1] = df['y0'].plot.hist(ax=axs[0,1],bins=50)
axs[0,1].axvline( x = qNom[0], ymin=0,ymax=20,color='red',linestyle='dashed')
axs[1,1] = df['y1'].plot.hist(ax=axs[1,1],bins=50)
axs[1,1].axvline( x = qNom[1], ymin=0,ymax=20,color='red',linestyle='dashed')
plt.show()



