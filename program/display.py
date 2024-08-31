from utils import *
from world import world
from data import *


class Image:
    def __init__(self):
        self.x = np.arange(dim.cell_size/2, x_size, dim.cell_size)
        self.y = np.arange(dim.cell_size/2, y_size, dim.cell_size)
        self.fig, self.ax = plt.subplots()

    def re_init(self):
        self.fig.clf()
        self.x = np.arange(dim.cell_size/2, x_size, dim.cell_size)
        self.y = np.arange(dim.cell_size/2, y_size, dim.cell_size)
        self.fig, self.ax = plt.subplots()

    @measure_execution_time
    def plot_function(self, values, plot_type):
        cmap = plt.colormaps["rainbow"]  # magma

        if plot_type == "binary":
            im = self.ax.pcolormesh(self.x, self.y, values, cmap=cmap, shading='nearest', vmin=B_MIN, vmax=B_MAX)
            self.fig.colorbar(im, ax=self.ax, orientation='horizontal', label="Débit Binaire [Gb/s]")
        else:
            im = self.ax.pcolormesh(self.x, self.y, values, cmap=cmap, shading='nearest', vmin=-90, vmax=-40)
            self.fig.colorbar(im, ax=self.ax, orientation='horizontal', label="Power [dbm]")

        im.set_mouseover(True)

        self.ax.set_title(f"Réception pour {n} émetteur{'s' if n > 1 else ''}")
        self.ax.set_xlabel("x[m]")
        self.ax.set_ylabel("y[m]")
        self.ax.set_aspect('equal')
        world.draw_walls(self.ax)

    def draw_rays(self, rays):
        rays.draw_rays_mpl(self.ax)

    def plot_emitters(self, txs):
        for i in range(int(len(txs)/2)):
            self.ax.scatter(txs[2*i], txs[2*i+1], s=20, color='blue', marker="+")

    @staticmethod
    def show():
        plt.show()

    def extract(self, filename):
        self.fig.savefig(filename, format='png')


image = Image()
