import pandas as pd
import xarray as xr
import geopandas
import seaborn as sns

import matplotlib.pyplot as plt


def loadNC(path):
    ds = xr.open_dataset(path)
    return ds.to_dataframe()




def map_plot(df, grid, filename="map.png"):
    df_grid = pd.concat([grid.reset_index(drop=True), df.reset_index(drop=True)], axis=1, ignore_index=False)

    #gdf = geopandas.GeoDataFrame(df_grid, geometry=geopandas.points_from_xy(df_grid.lon, df_grid.lat))
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    ax = world.plot(color='white', edgecolor='black')
    # gdf.plot(ax=ax, color='red')

    # ax = geoplot.kdeplot(ssh_grid.head(1000), shade=True, cmap='Reds', projection=geoplot.crs.AlbersEqualArea())
    # geoplot.polyplot(ssh_grid, ax=ax, zorder=1)

    print("KDE")
    start_index = 4_000_000
    print(df_grid.iloc[start_index:start_index + 3_000_000])
    sns.kdeplot(data=df_grid,  # .iloc[start_index:start_index+3_000_000],
                x='lon',
                y='lat',
                fill=True,
                cmap='coolwarm',
                alpha=0.3,
                gridsize=20,
                levels=50,
                ax=ax)

    plt.savefig(filename)
    plt.show()
ssh = loadNC('../map_summarized.nc')
unod = loadNC('../reading_u/map_summarized.nc')
mesh = xr.open_dataset("fesom.mesh.diag.nc")
grid = pd.concat([mesh["lon"].to_dataframe(), mesh["lat"].to_dataframe()], axis=1)


#map_plot(ssh, grid, filename="shh_map.png")
map_plot(unod, grid, filename="unod_map.png")
