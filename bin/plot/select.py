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
	
	df["index"] = df.Status
	df = df.set_index("index")
	n = df.loc["Unique", "Reads"] + df.loc["Included", "Reads"]
	df = df.reset_index()
	
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
			"size": 8,
			"color": "black"
		}
		ax.annotate(**args)
	
	ax.set_ylim(0, df.Reads.max() + df.Reads.max()/8)
	ax.tick_params(axis="x", which="major", labelsize=12)
	suptitle = "({:,} unambiguous/{:,} reads)".format(n, df.Reads.sum())
	ax.set_title( title + "\n" + suptitle )
	
	fig.tight_layout()

	return fig
	############################################################################

##########################
if __name__ == "__main__":

	csv_path = sys.argv[1]
	base_path = sys.argv[2]

	names = ["Process", "Sample", "Status", "Reads"]
	df = pd.read_csv(csv_path, header=None, names=names)
	
	names = {
		("UNIQUE", "SAME_GENE"): "Unique",
		("INCLUDED", "SAME_GENE"): "Included",
		("EXCLUDED", "SAME_GENE"): "Excluded\nSame gene",
		("EXCLUDED", "DIFFERENT_GENE"): "Excluded\n Different gene",
		("UNRESOLVED", "UNRESOLVED"): "Unresolved",
	}
	
	df = df\
		.join(
			df.Status\
				.str.split("/")\
				.apply(pd.Series)\
				.rename(columns={0:"CountStatus", 1:"GeneStatus"})
		)\
		.loc[:,["CountStatus", "GeneStatus", "Reads"]]\
		.set_index(["CountStatus", "GeneStatus"])\
		.reindex(names.keys(), fill_value=0)
	
	df.index = df.index.to_series().map(names)
	df.index.name = "Status"
	df = df.reset_index()
	df["Status"] = pd.Categorical(df.Status, categories=names.values())
	
	plt = count_plot(df, "Selecting multi-mapped UMIs (majority vote)")
	plt.savefig(f"{base_path}.png")
	plt.savefig(f"{base_path}.pdf")

