import pandas as pd
import xarray as xr
import geopandas
import seaborn as sns
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
    plt.tricontourf(df_grid.lon, df_grid.lat,df_grid.sea_surface_elevation, cmap='coolwarm', alpha = 0.3)
    plt.colorbar()
    #plt.scatter(x=df_grid.lon, y=df_grid.lat, c = df_grid.sea_surface_elevation)
    # # sns.kdeplot(data=df_grid.head(10000),  # .iloc[start_index:start_index+3_000_000],
    #             x='lon',
    #             y='lat',
    #             fill=True,
    #             cmap='coolwarm',
    #             alpha=0.3,
    #             gridsize=50,
    #             levels=20,
    #             cbar=True,
    #             hue="sea_surface_elevation",
    #             warn_singular=False,
    #             ax=ax)

    #plt.savefig(filename, dpi=1000)
    plt.savefig("prova.png")
    plt.show()


ssh = loadNC('../map_summarized.nc')
unod = loadNC('../reading_u/map_summarized.nc')
vnod = loadNC("../reading_v/map_summarized_vnod.nc")
mesh = xr.open_dataset("fesom.mesh.diag.nc")
grid = pd.concat([mesh["lon"].to_dataframe(), mesh["lat"].to_dataframe()], axis=1)

print(ssh.head(10_000))
#ssh.head(1_000_000).plot()
#plt.show()
# print(unod.head(1000))
# unod.head(1_000_000).plot()
# plt.show()
# print(vnod.head(10))
map_plot(ssh, grid, filename="shh_map.png")
# map_plot(unod, grid, filename="unod_map.png")
# map_plot(vnod, grid, filename="vnod_map.png")
