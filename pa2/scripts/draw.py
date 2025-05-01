import numpy as np
import matplotlib.pyplot as plt
import argparse

def argparse_args():
    parser = argparse.ArgumentParser(description='Draw a floorplan.')
    parser.add_argument('-t', '--title', type=str, default='Floorplan', help='Title of the plot')
    parser.add_argument('-b', '--block', type=str, required=True, help='Path to the block file')
    parser.add_argument('-r', '--rpt', type=str, required=True, help='Path to the rpt file')
    parser.add_argument('-o', '--output', type=str, default='floorplan.png', help='Output file name')
    return parser.parse_args()

def draw_block(name, x, y, w, h):
    pass

def draw_floorplan(block_file, rpt_file, output_file, title):
    # Read block file
    outline_width = 0
    outline_height = 0
    n_blocks = 0
    n_terminals = 0
    blocks = []
    terminals = []
    with open(block_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            if parts[0] == 'Outline:':
                outline_width = int(parts[1])
                outline_height = int(parts[2])
            elif parts[0] == 'NumBlocks:':
                n_blocks = int(parts[1])
            elif parts[0] == 'NumTerminals:':
                n_terminals = int(parts[1])
            elif len(parts) == 3 and 'terminal' not in parts:
                # This is a block definition: name width height
                block_name = parts[0]
                block_width = int(parts[1])
                block_height = int(parts[2])
                blocks.append({
                    'name': block_name,
                    'width': block_width,
                    'height': block_height
                })
            elif 'terminal' in parts:
                # This is a terminal definition: name terminal x y
                terminal_name = parts[0]
                terminal_x = int(parts[2])
                terminal_y = int(parts[3])
                terminals.append({
                    'name': terminal_name,
                    'x': terminal_x,
                    'y': terminal_y
                })

    # Read rpt file
    chip_width = 0
    chip_height = 0
    blocks_with_coords = []
    with open(rpt_file, 'r') as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 3:
                # First line is the chip size
                parts = line.split()
                chip_width = int(parts[0])
                chip_height = int(parts[1])
            elif i > 4:
                # Subsequent lines are block coordinates
                parts = line.split()
                block_name = parts[0]
                x1 = int(parts[1])
                y1 = int(parts[2])
                x2 = int(parts[3])
                y2 = int(parts[4])
                blocks_with_coords.append({
                    'name': block_name,
                    'x1': x1,
                    'y1': y1,
                    'x2': x2,
                    'y2': y2,
                })
    
    # Draw the blocks and terminals
    plt.figure(figsize=(10, 10))
    plt.xlim(0, outline_width)
    plt.ylim(0, outline_height)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(title)
    plt.plot([0, outline_width], [0, 0], 'k-', linewidth=1)
    plt.plot([0, 0], [0, outline_height], 'k-', linewidth=1)
    plt.plot([outline_width, outline_width], [0, outline_height], 'k-', linewidth=1)
    plt.plot([0, outline_width], [outline_height, outline_height], 'k-', linewidth=1)
    # draw chip outline with light blue
    plt.plot([chip_width, chip_width], [0, chip_height], 'b-', linewidth=2)
    plt.plot([0, chip_width], [chip_height, chip_height], 'b-', linewidth=2)
    for block in blocks_with_coords:
        x1 = block['x1']
        y1 = block['y1']
        x2 = block['x2']
        y2 = block['y2']
        plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], 'k-', linewidth=1)
        plt.fill_between([x1, x2], y1, y2, color='gray', alpha=0.5)
        plt.text((x1 + x2) / 2, (y1 + y2) / 2, block['name'], horizontalalignment='center', verticalalignment='center')
    plt.savefig(output_file, dpi=100)
    plt.close()
        

if __name__ == "__main__":
    args = argparse_args()
    draw_floorplan(args.block, args.rpt, args.output, args.title)