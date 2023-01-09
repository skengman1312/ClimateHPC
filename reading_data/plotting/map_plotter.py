import pandas as pd
import xarray as xr
import geopandas
import seaborn as sns

import matplotlib.pyplot as plt


def loadNC(path):
    ds = xr.open_dataset(path)
    return ds.to_dataframe()


ssh = loadNC('../map_summarized.nc')
unod = loadNC('../reading_u/map_summarized.nc')
mesh = xr.open_dataset("fesom.mesh.diag.nc")
grid = pd.concat([mesh["lon"].to_dataframe(), mesh["lat"].to_dataframe()], axis=1)

ssh_grid = pd.concat([grid.reset_index(drop=True), ssh.reset_index(drop=True)], axis=1, ignore_index=False)

gdf = geopandas.GeoDataFrame(ssh_grid, geometry=geopandas.points_from_xy(ssh_grid.lon, ssh_grid.lat))
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
ax = world.plot(color='white', edgecolor='black')
# gdf.plot(ax=ax, color='red')

# ax = geoplot.kdeplot(ssh_grid.head(1000), shade=True, cmap='Reds', projection=geoplot.crs.AlbersEqualArea())
# geoplot.polyplot(ssh_grid, ax=ax, zorder=1)

print("KDE")
sns.kdeplot(data=ssh_grid.head(3000000),
            x='lon',
            y='lat',
            fill=True,
            cmap='coolwarm',
            alpha=0.3,
            gridsize=200,
            levels=20,
            ax=ax)

plt.show()
print(ssh)
print(grid)
print(ssh_grid)
print(gdf.head())
