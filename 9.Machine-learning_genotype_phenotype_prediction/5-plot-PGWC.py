%matplotlib inline
import os, pickle, string
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats

from statannotations.Annotator import Annotator

font_path = "/home/zhoujb/local/font/Times New Roman.ttf"
mpl.font_manager.fontManager.addfont(font_path)
prop = mpl.font_manager.FontProperties(fname=font_path)
mpl.rcParams['font.family'] = prop.get_name()
mpl.rcParams['font.sans-serif'] = ["Times New Roman"]
mpl.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42

TEST_RES_PATH = "/data2/zhoujb/project/hhf/250718/Result/testRes4gene2K/"
PLOT_RES_PATH = "/data2/zhoujb/project/hhf/250718/plot4gene2K/"

model_list = ["SVR",  "EN", "Ridge", "BR", "ET"]
pre_list = ["10w", "150gene", "4gene2k"]

data_plot = pd.DataFrame()
for model_name in model_list:
    for pre_name in pre_list:
        with open(os.path.join(TEST_RES_PATH, "{}-{}-PGWC_BLUE.pickle".format(model_name, pre_name)), "rb") as in_f:
            tmp_data = pickle.load(in_f)
            tmp_data = tmp_data[tmp_data["Type"]=="R"].copy()
            tmp_data = tmp_data.groupby(["Model", "xName"])["Score"].mean().to_frame().reset_index()
            data_plot = pd.concat([data_plot, tmp_data], axis=0, ignore_index=True)

data_plot.head()

sns.set(font_scale=1.5, style="ticks")
plt.figure(figsize=(8,6))

markers = ['D', 'o', 'v']
predictor_name = data_plot["xName"].unique()
marker_dict = {k:v for k, v in zip(predictor_name, markers)}

ax = sns.boxplot(x="Model", y="Score", data=data_plot, color="white", showfliers=False)
# Add transparency to colors
for patch in ax.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .3))
    
#my_palette = sns.color_palette("husl", 3)
my_palette = ["#F73858", "#008080", "#F1D6D4"]
color_num = 0
for index in data_plot.index:
    if color_num == 3:
        color_num = 0
        
    plt.plot(data_plot.loc[index, "Model"], 
             data_plot.loc[index, "Score"], 
             marker_dict[data_plot.loc[index, "xName"]], 
             label=data_plot.loc[index, "xName"],
             fillstyle='none', color=my_palette[color_num],
             markersize=16, mew=3)
    
    color_num += 1
    
plt.ylabel('Pearson\'s R of PGWC')
plt.xlabel(None)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles[-3:], labels=labels[-6:], loc='upper right', 
          bbox_to_anchor=(1.3, 1), frameon=False, prop={'size': 16})

plt.savefig("../PGWC_BLUE.pdf", format="pdf", bbox_inches='tight', transparent=True)

