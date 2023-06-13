#Displays the up and down regulated differentially expressed genes (DEGS) per dataset as a bar plot.

d1_up = 4194 
d1_down = 2723 

d3_up = 1731
d3_down = 1936

d3_stage_up = 258
d3_stage_down = 717

d3_subtype_up = 432
d3_subtype_down = 571

import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 15})

plotdata = pd.DataFrame({

    "Up-regulated":[4194,0,1731,258,432],

    "Down-regulated":[2723,0,1936,717,571]},

    index=["D1", "D2", "D3", "D3_stage", "D3_subtype"])

plotdata.plot(kind="bar",figsize=(15, 8), color=['teal', 'tomato'])

#plt.title("Number of DEGs per dataset")

#plt.grid()

plt.xlabel("Datasets")

plt.ylabel("Number of DEGs")


#Venn diagram to show the intersection between datasets, given the values

#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
%matplotlib inline

venn2(subsets = (3607, 1144, 587), set_labels = ('Dataset 1', 'Dataset3'), set_colors=('red', 'teal'), alpha = 0.5)
venn2_circles(subsets = (3607, 1144, 587), linestyle='dashed', linewidth= 1 , color='k')
