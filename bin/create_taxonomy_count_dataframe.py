#!/usr/bin/env python3
__doc__ = """Parse Kaiju output to produce mapped read count dataframe."""
__author__ = "Fredrik Boulund"
__year__ = "2017"

from sys import argv, exit
import os
import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=". ".join([__doc__, __author__, __year__]))

    parser.add_argument("KAIJU", nargs="+",
            help="Kaiju summary output files.")
    parser.add_argument("-l", "--lineage", dest="lineage",
            default="superkingdom,phylum,class,order,family,genus,species",
            help="Extract lineage from Kaiju output (produced with -l). [%(default)s].")
    parser.add_argument("-o", "--output", 
            default="kaiju_dataframe.csv",
            help="Output filename [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_kaiju_summary_file(filename):
    """Read Kaiju summary file into Pandas DataFrame.

    Assumes filename suffix is the taxonomic level, e.g.
    samplename_visit.kaiju.summary.genus
    """
    sample_name = filename.split(".")[0]
    level = filename.split(".")[-1]
    df = pd.read_table(filename,
            index_col=2,
            skiprows=2,
            skipfooter=5,
            engine="python", # Required to use skipfooter
            dtype={sample_name: int},
            names=["Percent", sample_name, level])
    return df[sample_name]


def main(kaiju_files, output_filename, lineage):
    sample_dfs = []
    for kaiju_file in kaiju_files:
        sample_dfs.append(parse_kaiju_summary_file(kaiju_file))

    combined_df = pd.concat(sample_dfs, axis=1)

    ranks = lineage.split(",")
    if len(ranks) > len(combined_df.index[0].split(";")) - 1:
        print("WARNING: the lineage argument (-l) contains more ranks than the input data:\n", 
                lineage, "\n", 
                combined_df.index[0])
    split_ranks = (r.split(";")[:len(ranks)] for r in combined_df.index)
    rank_indices = list(zip(*split_ranks))
    lineage_index = pd.MultiIndex.from_arrays(rank_indices, names=ranks)
    combined_df.set_index(lineage_index, inplace=True)
    combined_df.to_csv(output_filename)
    return combined_df


if __name__ == "__main__":
    options = parse_args()
    main(options.KAIJU, options.output, options.lineage)
