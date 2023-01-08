import pandas as pd
import xarray as xr


def loadNC(path):
    ds = xr.open_dataset(path)
    return ds.to_dataframe()


print(loadNC('../map_summarized.nc'))
print(loadNC('../reading_u/map_summarized.nc'))

ds = xr.open_dataset('../reading_u/map_summarized.nc')
df = ds.to_dataframe()
print(df)
