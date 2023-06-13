d1_up = 4194 #5039
d1_down = 2723 #3541

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
