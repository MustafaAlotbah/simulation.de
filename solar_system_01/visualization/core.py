"""
    Core of the visualization module for simulation.de visualizations
    Mustafa Alotbah 2022
"""

import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt

from solar_system_01.simulation import Simulation


# Visualization class
class Visualization:
    def __init__(self, simulation: Simulation, bounds=None, simulation_steps=2501, foreground_color=None, background_color=None):
        self.foreground_color = foreground_color
        self.background_color = background_color

        if self.foreground_color is None:
            self.foreground_color = mpl.rcParams["text.color"]
        if self.background_color is None:
            self.background_color = mpl.rcParams["axes.labelcolor"]

        self.grid_color = '#181818'

        mpl.rcParams["text.color"] = self.foreground_color
        mpl.rcParams["axes.labelcolor"] = self.background_color
        mpl.rcParams["xtick.color"] = self.foreground_color
        mpl.rcParams["ytick.color"] = self.foreground_color

        self.bounds = bounds
        if bounds is None:
            self.bounds = [-5, 5, -10, 10]

        # Simulation to be visualized (or animated)
        self.simulation = simulation

        self.fig, self.ax1 = plt.subplots(1, 1, tight_layout=True)
        self.fig.set_facecolor(self.background_color)
        self.fig.set_size_inches(9, 16, forward=True)

        self.ax1.set_title("Space")
        self.ax1.set_facecolor(self.background_color)

        self.ax1.legend(
            handles=self.simulation.get_legend_handles(),
            framealpha=0.3
        )

        self.ax1.axis(self.bounds)
        self.ax1.set_aspect(1)

        self.ax1.set_xticks(range(self.bounds[0], self.bounds[1]))
        self.ax1.set_yticks(range(self.bounds[2], self.bounds[3]))
        self.ax1.grid(c=self.grid_color, zorder=-10)

        for artist in self.simulation.get_artists():
            self.ax1.add_patch(artist)

        def anim_update(_):
            simulation.update_artists()
            simulation.simulation_step()
            return simulation.get_artists()

        self.animation = mpl.animation.FuncAnimation(
            fig=self.fig,
            func=anim_update,
            interval=10,
            blit=True,
            frames=simulation_steps,
            repeat=False
        )

    @staticmethod
    def show():
        plt.show()

    def save_video(self, file_name):
        self.animation.save(file_name, fps=60, extra_args=['-vcodec', 'libx264'],
                 savefig_kwargs={'facecolor': self.background_color})

