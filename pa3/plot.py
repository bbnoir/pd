import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

def read_boundary_data(filename):
    """Read boundary coordinates from file"""
    coords = np.loadtxt(filename)
    return coords

def read_module_data(filename):
    """Read module coordinates from file with blank line separators"""
    with open(filename, 'r') as f:
        content = f.read().strip()
    
    modules = []
    for block in content.split('\n\n'):
        if block.strip():
            coords = np.array([list(map(float, line.split())) for line in block.strip().split('\n')])
            modules.append(coords)
    
    return modules

def read_wirelength(filename):
    """Read wirelength from info file"""
    with open(filename, 'r') as f:
        wirelength = float(f.readline().strip())
    return wirelength

def main():
    if len(sys.argv) != 4:
        print("Usage: python plot.py boundary_file modules_file info_file")
        sys.exit(1)
    
    boundary_file = sys.argv[1]
    modules_file = sys.argv[2]
    info_file = sys.argv[3]
    
    # Read data
    boundary = read_boundary_data(boundary_file)
    modules = read_module_data(modules_file)
    wirelength = read_wirelength(info_file)
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot boundary
    ax.plot(boundary[:, 0], boundary[:, 1], 'b-', linewidth=2)
    
    # Plot modules
    module_collection = PolyCollection(modules, facecolor='none', edgecolor='red', linewidth=1)
    ax.add_collection(module_collection)
    
    # Set plot properties
    ax.set_aspect('equal')
    ax.set_title(f'Placement Result (wirelength = {wirelength:.2f})')
    ax.autoscale_view()
    ax.grid(True)
    
    plt.tight_layout()
    # plt.show()
    plt.savefig('placement_result.png', dpi=300)

if __name__ == "__main__":
    main()
