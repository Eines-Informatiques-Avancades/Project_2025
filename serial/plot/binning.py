#!/usr/bin/python3

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


# Parse arguments.
def parse_args():
    parser = argparse.ArgumentParser(
        description = "Plot Binning from a data file."
    )
    parser.add_argument(
        'filename',
        type = str,
        help = "Path to the data file with the Binning data."
    )
    return parser.parse_args()

# Plot the data from the data file.
def plot_binning(input_file):
    # Read input file, ignoring lines starting with '#'.
    df = pd.read_csv (
        input_file,
        comment = '#',
        sep='\s+',
        header = None,
        names = ['ColName', 'BinSize', 'Mean', 'Variance', 'StdDev']
    )

    # Group data in terms of the first column (Measure type).
    groups = df.groupby('ColName')
    n_groups = len(groups)

    # Create a single file, with a subplot for each measure.
    fig, axs = plt.subplots(
        n_groups, 1, figsize=(8, 4 * n_groups), sharex=True
    )

    # In case there's only one measure
    if n_groups == 1:
        axs = [axs]

    for ax, (col_name, group) in zip(axs, groups):
        ax.plot(group['BinSize'], group['StdDev'], marker='o', linestyle='-')
        ax.set_title(f'{col_name}')
        ax.set_xlabel('  ')
        ax.set_ylabel('Standard Deviation')
        ax.grid(True)

    output_file = os.path.splitext(input_file)[0] + '.pdf'
    plt.savefig(output_file, format = 'pdf')
    plt.close()

    print(f"Binning plot saved as {output_file}")


if __name__ == '__main__':
    args = parse_args()
    plot_binning(args.filename)
