import numpy as np
import matplotlib.pyplot as plt

ms_3D_no_glyphs_10 = np.array([3.74, 4.69, 11.75, 63.97])
ms_2D_no_glyphs_10 = np.array([2.58, 3.26, 5.67, 26.23])
segment_count_10_traj = np.array([100,1000,10000,100000])
traj_count_10 = np.array([10,100,1000,10000])


ms_3D_no_glyphs_100 = np.array([3.57,5.45, 13.90, 95.94])
ms_2D_no_glyphs_100 = np.array([ 2.70,3.60, 8.06, 49.12])
segment_count_100_traj = np.array([1000,10000,100000,1000000])
traj_count_100 = np.array([10,100,1000,10000])


ms_3D_no_glyphs_1000 = np.array([4.68,12.09,87.84])
ms_2D_no_glyphs_1000 = np.array([3.14,6.65,44.08])
segment_count_1000_traj = np.array([10000,100000,1000000])
traj_count_1000 = np.array([10,100,1000])


# Helper function to plot with variable marker size
def plot_with_count(ax, x, y, traj_counts, color, marker, linestyle, label):    
    ax.plot(x, y, color=color, marker=marker, linestyle=linestyle, label=label, alpha=0.7)    
    # Annotate each point with its trajectory count
    for i in range(len(x)):
        ax.text(x[i], y[i], f"{traj_counts[i]}", fontsize=9, ha='right', va='bottom')


# Create a figure with 1 row, 2 columns
#fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
fig, ax = plt.subplots(figsize=(8, 6))
# Plot 3D
plot_with_count(ax, segment_count_10_traj, ms_2D_no_glyphs_10, traj_count_10, 'r','o',':', label="2D, Nodes=10")
plot_with_count(ax, segment_count_100_traj, ms_2D_no_glyphs_100, traj_count_100, 'g', 's',':', label="2D, Nodes=100")
plot_with_count(ax, segment_count_1000_traj, ms_2D_no_glyphs_1000, traj_count_1000, 'b','d',':', label="2D, Nodes=1000")


# Plot 2D
plot_with_count(ax, segment_count_10_traj, ms_3D_no_glyphs_10, traj_count_10, 'r','o','-', label="3D, Nodes=10")
plot_with_count(ax, segment_count_100_traj, ms_3D_no_glyphs_100, traj_count_100, 'g', 's','-', label="3D, Nodes=100")
plot_with_count(ax, segment_count_1000_traj, ms_3D_no_glyphs_1000, traj_count_1000, 'b','d', '-', label="3D, Nodes=1000")
ax.set_title("Frontside Tube Pipeline - Back&Frontside Tube Pipeline")
ax.set_xlabel("Segment Count")
ax.set_ylabel("Rendering Time (ms)")
ax.set_xscale("log")
ax.set_yscale("linear")
ax.legend()
ax.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()
