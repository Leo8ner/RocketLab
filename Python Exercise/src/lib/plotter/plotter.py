import numpy as np
import matplotlib.pyplot as plt

class ResultPlotter:
    """Class to plot a trajectory."""
    def __init__(self, results, title="Optimal Trajectory", gim_dyn=False, save=False):
        """
        Initialize the ResultPlotter with necessary data.

        Parameters:
        results (dict): Dictionary containing the results of the optimization.
        title (str): Title for the plot.
        """
        dt               = results['dt']  # Time steps
        X                = results['X']   # State trajectory
        X_dot            = results['X_dot']  # State derivatives
        U                = results['U']   # Control trajectory
        n                = dt.shape[0]  # Number of time steps
        self.time        = np.insert(np.cumsum(dt), 0, 0.0)  # insert t=0 at start
        self.x           = X[0, :].reshape((n+1,-1))  # Horizontal position
        self.z           = X[1, :].reshape((n+1,-1))  # Vertical position
        self.theta       = X[2, :].reshape((n+1,-1))  # Pitch angle
        self.vx          = X[3, :].reshape((n+1,-1))  # Horizontal velocity
        self.vz          = X[4, :].reshape((n+1,-1))  # Vertical velocity
        self.wy          = X[5, :].reshape((n+1,-1))  # Angular velocity
        self.vx_dot      = X_dot[3, :].reshape((n,-1))
        self.vz_dot      = X_dot[4, :].reshape((n,-1))
        self.wy_dot      = X_dot[5, :].reshape((n,-1))
        self.T           = U[0, :].reshape((n,-1))  # Engine thrust
        self.gimbal      = U[1, :].reshape((n,-1))  # Gimbal angle
        self.gim_real    = None  # Real gimbal angle, if dynamic gimbal is used
        if gim_dyn:
            self.gim_real = X[6, :].reshape((n+1,-1))
        self.title       = title
        self.save         = save  # Flag to save the plots
        self._plot_results()  

    def _plot_results(self):
        """
        Generate and display the plots for the given data.
        """
        fig, axes = plt.subplots(2, 4, figsize=(24, 10))
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

        # Plot linear accelerations
        axes[2].axhline(0, color='black', linestyle='--')
        axes[2].set_xlabel("Time (s)")
        axes[2].set_ylabel("Acceleration (m/s²)")
        axes[2].set_title("Acceleration vs Time")
        axes[2].plot(self.time[:-1], self.vx_dot, label="ax", color='orange')
        axes[2].plot(self.time[:-1], self.vz_dot, label="az", color='purple')
        axes[2].legend()

        # Plot pitch angle
        axes[3].axhline(0, color='black', linestyle='--')
        axes[3].set_xlabel("Time (s)")
        axes[3].set_ylabel("Angle (deg)")
        axes[3].set_title("Pitch Angle vs Time")
        axes[3].plot(self.time, self.theta * 180 / np.pi, label="theta", color='orange')
        axes[3].legend()

        # Plot angular velocity
        axes[4].axhline(0, color='black', linestyle='--')
        axes[4].set_xlabel("Time (s)")
        axes[4].set_ylabel("Angular rate (deg/s)")
        axes[4].set_title("Pitch Rate vs Time")
        axes[4].plot(self.time, self.wy * 180 / np.pi, label="wy", color='purple')
        axes[4].legend()

        # Plot angular acceleration
        axes[5].axhline(0, color='black', linestyle='--')
        axes[5].set_xlabel("Time (s)")
        axes[5].set_ylabel("Angular acceleration (deg/s²)")
        axes[5].set_title("Pitch Acceleration vs Time")
        axes[5].plot(self.time[:-1], self.wy_dot * 180 / np.pi, label="wy_dot", color='orange')
        axes[5].legend()

        # Plot engine thrust commands
        axes[6].axhline(0, color='black', linestyle='--')
        axes[6].set_xlabel("Time (s)")
        axes[6].set_ylabel("Thrust (N)")
        axes[6].set_title("Engine's Thrust vs Time")
        axes[6].plot(self.time[:-1], self.T, label="T", color='orange')
        axes[6].legend()

        # Plot gimbal angle
        axes[7].axhline(0, color='black', linestyle='--')
        axes[7].set_xlabel("Time (s)")
        axes[7].set_ylabel("Angle (deg)")
        axes[7].set_title("Gimbal Angle vs Time")
        if self.gim_real is not None:
            axes[7].plot(self.time[:-1], self.gimbal * 180 / np.pi, label="gimbal_cmd", color='purple', linestyle='--')
            axes[7].plot(self.time, self.gim_real * 180 / np.pi, label="gimbal", color='orange')
        else:
            axes[7].plot(self.time[:-1], self.gimbal * 180 / np.pi, label="gimbal", color='purple')
        axes[7].legend()
        # Save plots if required
        if self.save:
            plt.savefig(f"{self.title.replace(' ', '_').lower()}.png", dpi=300, bbox_inches='tight')
            print(f"Plot saved as {self.title.replace(' ', '_').lower()}.png")
            
        # Adjust layout for better spacing
        fig.tight_layout()
        plt.show()
