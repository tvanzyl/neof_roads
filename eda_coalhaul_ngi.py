import pandas as pd
import numpy as np
import os
import json
import imageio
import rasterio
import geopandas

from rasterio import features
from rasterio import windows
from rasterio import transform

data_pth_ngi = '/media/tvanzyl/data/ngi/'
data_pth_csv = '/media/tvanzyl/data/neof/{filename}'

gdf_18 = geopandas.read_file(data_pth_csv.format(filename='coalhaul_visuals_2018.csv'))
gdf_18[gdf_18==''] = None
gdf_18 = gdf_18.set_geometry( geopandas.GeoSeries.from_wkt(gdf_18.geom).set_crs("EPSG:4326") )
gdf_18.assessment_date = pd.to_datetime(gdf_18.assessment_date)
gdf_18.drop_duplicates(['road_nr', 'link_nr', 'from_km'], inplace=True, ignore_index=True)
gdf_18.set_index(['road_nr', 'link_nr', 'from_km',], inplace=True)

cols = ['geometry', 
        'risfsa_class', 'moisture_zone', 'terrain', 'functional_class', 'road_type', 'road_type_coto', 
        'road_width_m', 
        'width_m', 'lanes', 
        'surface_failures_d', 'surface_failures_e', 
        'surface_cracks_d', 'surface_cracks_e', 
        'stone_loss_d', 'stone_loss_e', 
        # 'dry_brittle_d', 'dry_brittle_e', 
        'bleeding_d', 'bleeding_e', 
        'block_cracks_d', 'block_cracks_e', 'longitudinal_cracks_d', 'longitudinal_cracks_e', 'transverse_cracks_d', 'transverse_cracks_e', 'crocodile_creacks_d', 'crocodile_cracks_e', 
        'pumping_d', 'pumping_e', 
        'rutting_d', 'rutting_e', 
        'undulation_settlement_d', 'undulation_settlement_e', 
        'patching_d', 'patching_e', 
        'failures_potholes_d', 'failures_potholes_e', 
        'edge_breaks_d', 'edge_breaks_e', 
        'riding_quality',  
        'vci']

# 'assessment_date''route''stone_loss_an','shoulders_unpaved','surface_drainage',

gdf_18 = gdf_18[cols]

gdf_18.dropna(axis='columns', how='all', inplace=True)
gdf_18.dropna(axis='index', how='any', inplace=True)

index = gdf_18.index
geometry = gdf_18.loc[index].geometry

numeric_cols = ['surface_failures_d', 'surface_failures_e',
                'failures_potholes_d', 'failures_potholes_e', 
                # 'surface_cracks_d', 'surface_cracks_e', 
                'stone_loss_d', 'stone_loss_e',
                # 'dry_brittle_d', 'dry_brittle_e', 
                'bleeding_d', 'bleeding_e',
                'block_cracks_d', 'block_cracks_e', 'longitudinal_cracks_d', 'longitudinal_cracks_e', 
                'transverse_cracks_d', 'transverse_cracks_e', 'crocodile_creacks_d', 'crocodile_cracks_e', 
                'pumping_d', 'pumping_e', 'rutting_d', 'rutting_e', 
                'undulation_settlement_d', 'undulation_settlement_e', 
                'edge_breaks_d', 'edge_breaks_e', 
                'patching_d', 'patching_e', 
                'road_width_m','width_m', 'lanes', 'riding_quality', 'vci']

gdf_18 = gdf_18[numeric_cols].astype(float)
gdf = gdf_18.loc[index, numeric_cols]

tile_names = [name for name in os.listdir(data_pth_ngi) if os.path.isdir(os.path.join(data_pth_ngi, name))]

tile_name = tile_names[0]
tif_files = [f for f in os.listdir(data_pth_ngi+f"{tile_name}/") if f.endswith(".tif")]
file_name = tif_files[0]

src = rasterio.open(data_pth_ngi+f"{tile_name}/{file_name}")
gdf = gdf.set_geometry(geometry.to_crs(src.crs))
gdf = gdf.reset_index(drop=True)

i = 0
seq = []
idxs = []
for tile_name in tile_names:
    rstpath = data_pth_csv.format(filename=f'{tile_name}')

    try:
        os.mkdir(rstpath) 
        print(f"Directory '{rstpath}' created successfully.")
    except FileExistsError:
        print(f"Directory '{rstpath}' already exists.")    
    except FileNotFoundError:
        print(f"Parent directory not found for '{rstpath}'.")
    except PermissionError:
        print(f"You don't have permission to create '{rstpath}'.")

    tif_files = [f for f in os.listdir(data_pth_ngi+f"{tile_name}/") if f.endswith(".tif")]
    file_name = tif_files[0]

    for file_name in tif_files:
        src = rasterio.open(data_pth_ngi+f"{tile_name}/{file_name}")
        gdf_subset =  gdf.cx[src.bounds[0]:src.bounds[2],src.bounds[1]:src.bounds[3]].clip_by_rect(*src.bounds)

        for idx in gdf_subset.index:
            row = gdf_subset[idx]        
            bounds = row.bounds
            rst_rgb = src.read([1,2,3], window=windows.from_bounds(*bounds, src.transform))
            rst = rst_rgb[0]
            if rst.shape[0]*rst.shape[1] > 0:        
                imageio.v3.imwrite(data_pth_csv.format(filename=f'{tile_name}/{i}.tif'), rst_rgb)

                tmp_tran = transform.from_bounds(*bounds, rst.shape[1], rst.shape[0])
                tmp_msk = features.rasterize( ((row, 1),), fill=0, transform=tmp_tran, out_shape=rst.shape, all_touched=True)
                
                imageio.imwrite(data_pth_csv.format(filename=f'{tile_name}/{i}_mask.tif'), tmp_msk)

                tmp_seq = [rst[i[0],i[1]] for i in np.argwhere(tmp_msk)]
                seq.append((idx,tmp_seq))
                idxs.append(idx)
                
                i = i + 1

        print(f"Generated {i} Patches")

gdf.loc[idxs].reset_index().to_csv(data_pth_csv.format(filename=f'targets.csv'))

# with open(data_pth_csv.format(filename=f'{rstsource}/sequences.json'), 'w') as f:
#     json.dump(eval(str(seq)), f, ensure_ascii=False)