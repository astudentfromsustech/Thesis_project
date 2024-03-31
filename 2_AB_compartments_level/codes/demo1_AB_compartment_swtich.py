import pathlib as p
import pandas as pd

import numpy as np
import matplotlib

matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import math
import seaborn as sns
import sys
import networkx as nx
# Increase the recursion limit to handle deep recursions
sys.setrecursionlimit(10000)
# import typing as t
plt.style.use('ggplot')

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)

#### draw the pieplot
values = [7573, 42, 4008, 210]  # BB,BA,AA,AB

# colors = ['gray', '#B8860B', '#00008B', '#800000']
colors = ['gray','orange','green','lightblue']
fig, ax = plt.subplots(figsize=(5, 5), dpi=200)

# Removing the autopct parameter to avoid displaying percentages
plt.pie(values, colors=colors, startangle=90)
plt.axis('equal')

# Display the plot
plt.show()
