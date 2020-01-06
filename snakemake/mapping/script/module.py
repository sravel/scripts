#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @package module.py
# @author Sebastien Ravel

##################################################
## Modules
##################################################
## Python modules
import argparse, pysam, traceback, sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict, OrderedDict
import pandas as pd
from tempfile import NamedTemporaryFile

# environment settings:
pd.set_option('display.max_column',None)
pd.set_option('display.max_rows',None)
pd.set_option('display.max_seq_items',None)
pd.set_option('display.max_colwidth', 500)
pd.set_option('expand_frame_repr', True)
# pd.options.display.width = None


def parse_idxstats(files_list = None, out_csv = None, sep="\t"):
	from pathlib import Path
	from collections import defaultdict, OrderedDict
	import pandas as pd
	dico_mapping_stats = defaultdict(OrderedDict)
	for csv_file in files_list:
		sample = Path(csv_file).stem.split("_")[0]
		reference = Path(csv_file).parent.name
		print(reference)
		df = pd.read_csv(csv_file, sep="\t", header=None,names=["chr","chr_size","map_paired","map_single"], index_col=False)
		# print(df)
		unmap = df[df.chr == '*'].map_single.values[0]
		df = df[df.chr != '*']
		map_total = df["map_single"].sum()+df["map_paired"].sum()
		size_lib = map_total+unmap
		poucent = map_total/size_lib
		dico_mapping_stats[f"{reference} {sample}"]["size_lib"] = size_lib
		dico_mapping_stats[f"{reference} {sample}"]["map_total"] = map_total
		dico_mapping_stats[f"{reference} {sample}"]["poucent"] = f"{poucent*100:.2f}%"
	dataframe_mapping_stats = pd.DataFrame.from_dict(dico_mapping_stats, orient='index')
	with open(out_csv, "w") as out_csv_file:
		# print(f"Library size:\n{dataframe_mapping_stats}\n")
		dataframe_mapping_stats.to_csv(out_csv_file, index=True, sep=sep)
