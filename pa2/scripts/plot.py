#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse

def argparse_args():
    parser = argparse.ArgumentParser(description='Plot SA')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file')
    parser.add_argument('-o', '--output', type=str, default='plot.png', help='Output file')
    return parser.parse_args()

def plot_sa(input_file, output_file):
    # iter, cost, wl, area
    data = np.loadtxt(input_file, delimiter=',')

    iter = data[:, 0] # x
    cost = data[:, 1] # y1
    wl = data[:, 2] # y2
    area = data[:, 3] # y3
    temperature = data[:, 4] # y4

    # Create a figure with 4 subplots
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

    # Plot cost in the first subplot
    ax1.plot(iter, cost, 'b-')
    min_cost = np.min(cost)
    ax1.axhline(y=min_cost, color='blue', linestyle='--', alpha=0.7)
    ax1.set_ylabel('Cost', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.set_title('Cost over Iterations')
    ax1.grid(True)

    # Plot wirelength in the second subplot
    ax2.plot(iter, wl, 'r-')
    min_wl = np.min(wl)
    ax2.axhline(y=min_wl, color='red', linestyle='--', alpha=0.7)
    ax2.set_ylabel('Wirelength', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.set_title('Wirelength over Iterations')
    ax2.grid(True)

    # Plot area in the third subplot
    ax3.plot(iter, area, 'g-')
    min_area = np.min(area)
    ax3.axhline(y=min_area, color='green', linestyle='--', alpha=0.7)
    ax3.set_ylabel('Area', color='g')
    ax3.tick_params(axis='y', labelcolor='g')
    ax3.set_title('Area over Iterations')
    ax3.grid(True)
    
    # Plot temperature in the fourth subplot
    ax4.plot(iter, temperature, 'm-')
    ax4.set_xlabel('Iteration')
    ax4.set_ylabel('Temperature', color='m')
    ax4.tick_params(axis='y', labelcolor='m')
    ax4.set_title('Temperature over Iterations')
    ax4.grid(True)

    # Adjust layout
    plt.tight_layout()
    
    # Save the plot with higher resolution
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    

if __name__ == "__main__":
    args = argparse_args()
    plot_sa(args.input, args.output)