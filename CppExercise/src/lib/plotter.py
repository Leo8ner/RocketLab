import numpy as np
import matplotlib.pyplot as plt
import sys

class ResultPlotter:
    """Class to plot a trajectory."""
    def __init__(self, path, N_steps, N_X, title="Optimal Trajectory", gim_dyn=False, save=False):
        """
        Initialize the ResultPlotter with necessary data.

        Parameters:
        results (dict): Dictionary containing the results of the optimization.
        title (str): Title for the plot.
        """
        dt = np.loadtxt(path, delimiter=",", usecols=range(1), skiprows=1, max_rows=N_steps).T
        X  = np.loadtxt(path, delimiter=",", usecols=range(1, N_X+1), skiprows=1).T
        U  = np.loadtxt(path, delimiter=",", usecols=range(N_X+1, N_X+3), skiprows=1, max_rows=N_steps).T


        n                = N_steps  # Number of time steps
        self.time        = np.insert(np.cumsum(dt), 0, 0.0)  # insert t=0 at start
        self.x           = X[0, :].reshape((n+1,-1))  # Horizontal position
        self.z           = X[1, :].reshape((n+1,-1))  # Vertical position
        self.theta       = X[2, :].reshape((n+1,-1))  # Pitch angle
        self.vx          = X[3, :].reshape((n+1,-1))  # Horizontal velocity
        self.vz          = X[4, :].reshape((n+1,-1))  # Vertical velocity
        self.wy          = X[5, :].reshape((n+1,-1))  # Angular velocity
        self.T           = U[0, :].reshape((n,-1))  # Engine thrust
        self.gimbal      = U[1, :].reshape((n,-1))  # Gimbal angle
        self.gim_real    = None  # Real gimbal angle, if dynamic gimbal is used
        if gim_dyn:
            self.gim_real = X[6, :].reshape((n+1,-1))
        self.title       = title
        self.save        = save  # Flag to save the plots
        self._plot_results()  

    def _plot_results(self):
        """
        Generate and display the plots for the given data.
        """
        fig, axes = plt.subplots(2, 3, figsize=(24, 10))
        axes      = axes.flatten()

        plt.gcf().canvas.manager.set_window_title(self.title)

        # Plot positions
        axes[0].axhline(0, color='black', linestyle='--')
        axes[0].plot(self.time, self.x, color='orange', label="x")
        axes[0].plot(self.time, self.z, color='purple', label="z")
        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Position (m)")
        axes[0].set_title("Position vs Time")
        axes[0].legend()

        # Plot linear velocities
        axes[1].axhline(0, color='black', linestyle='--')
        axes[1].set_xlabel("Time (s)")
        axes[1].set_ylabel("Velocity (m/s)")
        axes[1].set_title("Velocity vs Time")
        axes[1].plot(self.time, self.vx, label="vx", color='orange')
        axes[1].plot(self.time, self.vz, label="vz", color='purple')
        axes[1].legend()


        # Plot pitch angle
        axes[2].axhline(0, color='black', linestyle='--')
        axes[2].set_xlabel("Time (s)")
        axes[2].set_ylabel("Angle (deg)")
        axes[2].set_title("Pitch Angle vs Time")
        axes[2].plot(self.time, self.theta * 180 / np.pi, label="theta", color='orange')
        axes[2].legend()

        # Plot angular velocity
        axes[3].axhline(0, color='black', linestyle='--')
        axes[3].set_xlabel("Time (s)")
        axes[3].set_ylabel("Angular rate (deg/s)")
        axes[3].set_title("Pitch Rate vs Time")
        axes[3].plot(self.time, self.wy * 180 / np.pi, label="wy", color='purple')
        axes[3].legend()

        # Plot engine thrust commands
        axes[4].axhline(0, color='black', linestyle='--')
        axes[4].set_xlabel("Time (s)")
        axes[4].set_ylabel("Thrust (N)")
        axes[4].set_title("Engine's Thrust vs Time")
        axes[4].plot(self.time[:-1], self.T, label="T", color='orange')
        axes[4].legend()

        # Plot gimbal angle
        axes[5].axhline(0, color='black', linestyle='--')
        axes[5].set_xlabel("Time (s)")
        axes[5].set_ylabel("Angle (deg)")
        axes[5].set_title("Gimbal Angle vs Time")
        if self.gim_real is not None:
            axes[5].plot(self.time[:-1], self.gimbal * 180 / np.pi, label="gimbal_cmd", color='purple', linestyle='--')
            axes[5].plot(self.time, self.gim_real * 180 / np.pi, label="gimbal", color='orange')
        else:
            axes[5].plot(self.time[:-1], self.gimbal * 180 / np.pi, label="gimbal", color='purple')
        axes[5].legend()
        # Save plots if required
        if self.save:
            plt.savefig(f"../output/{self.title.replace(' ', '_').lower()}.png", dpi=300, bbox_inches='tight')
            print(f"Plot saved as {self.title.replace(' ', '_').lower()}.png")
            
        # Adjust layout for better spacing
        fig.tight_layout()
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python plotter.py results.csv N_steps N_X [title] [gim_dyn] [save]")
        sys.exit(1)
    ResultPlotter("../output/" + sys.argv[1], int(sys.argv[2]), int(sys.argv[3]),
         title=sys.argv[4] if len(sys.argv) > 4 else "Optimal Trajectory",
         gim_dyn=bool(int(sys.argv[5])) if len(sys.argv) > 5 else False,
         save=bool(int(sys.argv[6])) if len(sys.argv) > 6 else False)
