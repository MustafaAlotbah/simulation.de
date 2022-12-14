import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches

#
class Particle:
    def __init__(self, x=0, y=0, vx=0, vy=0, mass=1, color='white', label="", is_fixed=False):
        # physical properties
        self.position = np.array([x, y], dtype=np.float)
        self.velocity = np.array([vx, vy], dtype=np.float)
        self.force = np.array([0.0, 0.0], dtype=np.float)
        self.mass = mass

        # other properties
        self.color = color
        self.label = label
        self.is_fixed = is_fixed
        self.apply_earth_gravity = False
        self.apply_gravitational_forces = True

        # what to draw on the plot
        self.artist = plt.Circle(
            tuple(self.position),
            radius=1,
            fc=self.color, ec=self.color,
            linewidth=5,
            zorder=10
        )

        # what to draw on the legend
        self.legend_handle = plt.Line2D(        # this is a -o- plot
            [0], [0],
            color='#ffffff00',          # this is the line's color (transparent)
            markerfacecolor=self.color,
            markeredgecolor=self.color,
            marker='o',
            markersize=10,
            label=self.label
        )

        # information to keep track of its position
        self.positions = [tuple(self.position)]

        _path = mpl_patches.Path(self.positions, [mpl_patches.Path.MOVETO])
        self.path_artist = mpl_patches.PathPatch(path=_path, ec=self.color, linewidth=0.75, fill=False, zorder=8, alpha=0.6)

    def get_legend_handle(self):
        return self.legend_handle

    def get_artists(self):
        return self.artist, self.path_artist

    def update_artist(self):
        # update the particle itself
        self.artist.center = self.position
        self.artist.set_color(self.color)
        self.artist.radius = max(1, np.sqrt(max(0, np.log10(self.mass))))/15     # size by design

        # update the particle's trace
        self.positions.append(tuple(self.position))
        self.positions = self.positions[-750:]          # keep the tail 750 points long at most
        _path = mpl_patches.Path(self.positions, [mpl_patches.Path.MOVETO] + [mpl_patches.Path.LINETO] * (len(self.positions) - 1))
        self.path_artist.set_path(_path)
