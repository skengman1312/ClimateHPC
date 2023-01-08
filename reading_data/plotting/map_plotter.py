import pandas as pd
import xarray as xr


def loadNC(path):
    ds = xr.open_dataset(path)
    return ds.to_dataframe()


print(loadNC('../map_summarized.nc'))
print(loadNC('../reading_u/map_summarized.nc'))
mesh = xr.open_dataset("fesom.mesh.diag.nc")
grid = pd.concat([mesh["lon"].to_dataframe(), mesh["lat"].to_dataframe()], axis=1)
print(grid)


