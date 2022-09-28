#!/usr/bin/python

import sys
from itertools import product
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib as mpl
mpl.rc('font', size=16)

###########################
def count_plot(df, title):#
###########################
	"""Plots counts. The argument data frame should have at least 2 columns
	named "Name" and "Count". The function returns a Figure object."""

	# annotation
	df["Percent"] = np.round(df.Reads / df.Reads.sum() * 100, 1)
	df["Annot"] = df.Reads.apply(lambda x: "{:,}".format(round(x,2)))
	df["Annot"] = df.Annot + " reads\n" + df.Percent.astype(str) + " %"

	# the plot
	fig = Figure(figsize=(8, 8))
	ax = fig.add_subplot(111)
	args = {
		"x": df.Status,
		"height": df.Reads,
		"color": "skyblue",
		"edgecolor": "black"
	}
	ax.bar(**args)

	# add number and percentage in the middle of the bar
	for i, annot in enumerate(df.Annot):
		args = {
			"s": annot,
			"xy": (i, df.iloc[i].Reads + df.Reads.max()/50),
			"ha": "center",
			"va": "bottom",
			"size": 16,
			"color": "black"
		}
		ax.annotate(**args)

	ax.set_ylim(0, df.Reads.max() + df.Reads.max()/8)
	ax.set_title(title)
	
	fig.tight_layout()

	return fig
	############################################################################

#csv_path = "results/temp_files/sample1/12_duplicates/sample1.pos_dups.csv"
#base_path = "tmp/test"

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]

	names = ["Process", "Sample", "Status", "Reads"]
	df = pd.read_csv(csv_path, header=None, names=names)
	df = df[["Status", "Reads"]]
	
	values = df.set_index("Status").Reads
	values.loc["Non duplicates"] = values.loc["Total"] - values.loc["Duplicates"]
	values = values.loc[["Non duplicates", "Duplicates"]].reset_index()
	status = ["Non duplicates", "Duplicates"]
	values["Status"] = pd.Categorical(values.Status, categories=status)
	
	plt = count_plot(values, "Position duplicates\n(same UMI and same mapping)")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

