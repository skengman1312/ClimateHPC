import pandas as pd
import xarray as xr
import geopandas
import seaborn as sns
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.animation as animation
from netCDF4 import Dataset, MFDataset
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
# import geoplot

import matplotlib.pyplot as plt


def loadNC(path):
    ds = xr.open_dataset(path)
    return ds.to_dataframe()


def map_plot(df, grid, filename="map.png"):
    df_grid = pd.concat([grid.reset_index(drop=True), df.reset_index(drop=True)], axis=1, ignore_index=False)

    # gdf = geopandas.GeoDataFrame(df_grid, geometry=geopandas.points_from_xy(df_grid.lon, df_grid.lat))
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    ax = world.plot(color='white', edgecolor='black')
    # # gdf.plot(ax=ax, color='red')
    #
    # ax = geoplot.kdeplot(df_grid.head(1000), shade=True, cmap='Reds', projection=geoplot.crs.AlbersEqualArea())
    # geoplot.polyplot(df_grid, ax=ax, zorder=1)
    # print(df_grid.sea_surface_elevation)
    print("KDE")

    # plt.tricontourf(df_grid.lon, df_grid.lat, df_grid.iloc[:, 2], cmap='coolwarm', alpha=0.3)
    plt.tricontourf(df_grid.lon, df_grid.lat, df_grid.iloc[:, 2], cmap='coolwarm', alpha=0.3, vmin=-1.90, vmax=2)

    # plt.scatter(x=df_grid.lon, y=df_grid.lat, c=df_grid.iloc[:,2], s=0.01, cmap='coolwarm')
    plt.colorbar()
    # # sns.kdeplot(data=df_grid.head(10000),  # .iloc[start_index:start_index+3_000_000],
    #             x='lon',
    #
    #             fill=True,
    #             cmap='coolwarm',
    #             alpha=0.3,
    #             gridsize=50,
    #             levels=20,
    #             cbar=True,
    #             hue="sea_surface_elevation",
    #             warn_singular=False,
    #             ax=ax)

    # sns.scatterplot(df_grid,  x='lon',  y='lat',  cmap='coolwarm',  hue="sea_surface_elevation",  ax=ax)
    plt.savefig(filename, dpi=1000)

    plt.show()


def animated_plot(df, grid):
    fig, ax = plt.subplots()
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

    def animate(i):
        ax.clear()
        world.plot(color='white', edgecolor='black', ax=ax)
        df_grid = pd.concat([grid.reset_index(drop=True), df.loc[[i]].reset_index(drop=True)], axis=1,
                            ignore_index=False)
        ax.tricontourf(df_grid.lon, df_grid.lat, df_grid.iloc[:, 2], cmap='coolwarm', alpha=0.3, vmin=-1.90, vmax=2)
        # plt.colorbar()
        plt.title('t = %i: ' % (i + 1))

    anim = animation.FuncAnimation(fig, animate, 12, interval=1000, blit=False)
    anim.save('us.gif', writer='imagemagick', fps=1, dpi=1000)
    # plt.show()


def round_plotting(i, kind="ssh", directory="final_plots", dpi = 1000):
    plt.ioff()
    d = {"ssh": {"title": "Sea Surface Elevation", "cbar": "elevation in meters", "fname": "map_ssh.png", "filepath":'../reading_ssh_v2/map_summarized.nc'},
         "unod": {"title": "Horizontal speed", "cbar": "speed in meters per second", "fname": "map_unod.png", "filepath":'map_summarized_unod.nc'},
         "vnod": {"title": "Vertical speed", "cbar": "speed in meters per second", "fname": "map_vnod.png", "filepath":'map_summarized_vnod.nc'}
         }
    mon = {0: "January", 1: "February", 2: "March", 3: "April", 4: "May", 5: "June",
                     6: "July", 7: "August", 8: "September", 9: "October", 10: "November", 11: "December"}

    f3 = Dataset('fesom.mesh.diag.nc')
    f5 = Dataset(d[kind]["filepath"])
    m = Basemap(projection='robin', lon_0=0, resolution='c')
    # print(f3)
    # print(f5.variables)
    v = list(f5.variables.keys())[0]

    # print(f"Variable: {v}")
    x, y = m(f3.variables['lon'], f3.variables['lat'])
    plt.figure(figsize=(15, 10))

    plt.tricontourf(x, y, f5.variables[v][i], cmap='coolwarm', extend='both', alpha=0.5, vmin=-1.90, vmax=2)
    cbar = plt.colorbar(orientation='horizontal', pad=0.03)
    cbar.set_label(d[kind]["cbar"])
    m.drawcoastlines()
    m.drawmapboundary(fill_color='0.9')
    m.fillcontinents()

    plt.title(f"{mon[i]} {d[kind]['title']}")
    plt.tight_layout()
    savename = f"{directory}/{mon[i].lower()}_{d[kind]['fname']}"
    plt.savefig(savename, dpi=dpi)


if __name__ == "__main__":
    ssh = loadNC('../reading_data/reading_ssh/map_summarized.nc')

    # ssh = loadNC('../map_summarized.nc')
    # unod = loadNC('../reading_u/map_summarized.nc')
    # vnod = loadNC("../reading_v/map_summarized_vnod.nc")
    # vnod.rename(columns={"speed": "vspeed"}, inplace=True)
    ssh.rename(columns={"sea_surface_elevation": "ssh"})
    mesh = xr.open_dataset("fesom.mesh.diag.nc")
    grid = pd.concat([mesh["lon"].to_dataframe(), mesh["lat"].to_dataframe()], axis=1)

    print("min: ", ssh.min())
    print("max ", ssh.max())

    # animated_plot(ssh, grid)

    #
    # cdata = pd.concat([unod, vnod], axis=1)       #ssh, unod, vnod], axis=1)
    #
    # # cdata.drop("time", axis=0)
    #
    # cdata.reset_index(1, inplace=True)
    # cdata.drop("time", axis=1, inplace=True)
    # book2 = KMeans(n_clusters=5, random_state=1312, verbose=True).fit_predict(cdata)
    # #counts = np.unique(book2.labels_, return_counts=True)
    # dbook = pd.DataFrame(book2)
    # print(dbook)
    # map_plot(dbook, grid, "cluster_map.png")
    # print(book2.shape)
    # print(counts)

    # ssh.head(1_000_000).plot()
    # plt.show()
    # print(unod.head(1000))
    # unod.head(1_000_000).plot()
    # plt.show()
    # print(vnod.head(10))
    # print(ssh.iloc[0])
    # map_plot(ssh.loc[[11]], grid, filename="shh_map.png")
    # map_plot(unod, grid, filename="unod_map.png")
    # map_plot(vnod, grid, filename="vnod_map.png")
    [[round_plotting(i, kind=t, dpi=200, directory="low_res_plots") for i in range(12)] for t in ["ssh","vnod","unod"]]

