# With help of ChatGPT, as I'm not an experienced Python programmer

import csv
import random
import numpy as np
import os

# Function to generate non-uniform random time points
def generate_time_points_non_uniform(count, start=0, end=60):
    # Using an exponential distribution for non-uniformity
    scale = (end - start) / 3  # Adjust scale for spread
    times = np.random.exponential(scale, count)
    times = start + (times / max(times)) * (end - start)  # Normalize to [start, end]
    
    return sorted(times)

# Function to generate random time points
def generate_time_points(count, start=0, end=60):
    return sorted(np.random.uniform(start, end, count))

# Function to fill positions based on time length
def generate_positions(length, axis):
    # Axis-specific tensor generation
    if axis == "x":
        return [[i/10.0, 0.0, 0.0] for i in range(length)]
    elif axis == "y":
        return [[0.0, i/10.0, 0.0] for i in range(length)]
    elif axis == "z":
        return [[0.0, 0.0, i/10.0] for i in range(length)]
    else:
        raise ValueError("Invalid axis. Choose 'x', 'y', or 'z'.")

def generate_tensors(length, axis):
    # Axis-specific tensor generation
    if axis == "x":
        return [[2.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5] for _ in range(length)]
    elif axis == "y":
        return [[0.5, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.5] for _ in range(length)]
    elif axis == "z":
        return [[0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 2.0] for _ in range(length)]
    else:
        raise ValueError("Invalid axis. Choose 'x', 'y', or 'z'.")

# Get the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define output file path in script's directory
output_file = os.path.join(script_dir, "DemoDebugTensors.csv")

# Header and data preparation
header = ["Timestamp", "ID", "Pos:0", "Pos:1", "Pos:2", "Diff:0", "Diff:1", "Diff:2", "Diff:3", "Diff:4", "Diff:5", "Diff:6", "Diff:7", "Diff:8"]
row_sum = 0
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)  # Write the header

    for iteration, axis in enumerate(["x", "y", "z"], start=0):
        # Number of random time points
        time_count = random.randint(100, 200)  # Random count between 100 and 200

        # Generate random time points
        random_times = generate_time_points(time_count)
        random_times.insert(0, 0)
        row_sum += len(random_times)
        current_traj_id = iteration

        # Generate positions for each time point
        random_positions = generate_positions(len(random_times), axis)

        # Generate tensors for each time point
        tensors = generate_tensors(len(random_times), axis)

        
        data = [
            [time, current_traj_id] + positions + tensors for time, positions, tensors in zip(random_times, random_positions, tensors)
        ]
        writer.writerows(data)  # Write the rows

print(f"CSV file '{output_file}' created successfully with {row_sum} rows.")