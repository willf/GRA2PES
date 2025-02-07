#!/usr/bin/env python
# coding: utf-8
# use xesmf_regrid_env
# Script to Regrid from CAMS domain to WRF Domain
# Colin Harkins - January 2020

# Use xesmf_regrid_env.yml to create python environment

import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
import ESMF
import pyproj
import os

# %% Inputs

# WRF file for WRF Domain
#out_example           = xr.open_dataset('wrfchemi_00z_d01.nc', chunks={'Time': 1,'emissions_zdim': 1})
input_dir             ='/wrk/charkins/emissions/CAMSv4.2/' 
out_example_fn        ='/wrk/d2/charkins/CAMSv4.2/template/template_00to12Z.nc'
out_example_fn2       ='/wrk/d2/charkins/CAMSv4.2/template/template_12to24Z.nc'
#out_dir               = '/wrk/d2/CAMSv4.2/regrid/'
out_dir               = '/wrk/users/charkins/emissions/CAMSv4.2/regrid/2019_fixedOM/Month01/'
out_filename_suffix   = ''
out_filename_type     = '.nc4'

os.chdir(input_dir)
# CAMS file for data
#CAMS_dir               = '/wrk/d2/CAMSv4.2/'
CAMS_dir               = '/wrk/d2/charkins/CAMSv4.2/'
CAMS_filename_prefix   = 'CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_'
CAMS_filename_suffix   ='_v4.2_monthly.nc'

# Species List from CAMS
CAMS_species = ['bc', 'ch4', 'co2', 'co', 'nh3', 'nmvoc', 'nox', 'oc', 'so2',
                 'voc1', 'voc2', 'voc3', 'voc4', 'voc5', 'voc6', 'voc7', 'voc8', 'voc9', 'voc10',
                 'voc11', 'voc12', 'voc13', 'voc14', 'voc15', 'voc16', 'voc17', 'voc18', 'voc19', 'voc20',
                 'voc21', 'voc22', 'voc23', 'voc24', 'voc25']
#CAMS_species = ['nox','ch4','co']
#WRF_species = ['NOX','CO','CO2']
WRF_species = ['NOX','CO','CO2','SO2','NH3','PM04','PM05','PM06','HC01',
		'HC02','HC05','HC06','HC07','HC08','HC09','HC10','HC11','HC14','HC18',
		'HC20','HC31','HC37','HC38','HC39','HC40','HC41','HC42','HC44','HC45',
                'HC48','HC49','HC50']

# Offroad and onroad (tnr, tro) need to be done separately due to different DOW profiles
#CAMS_sectors = ['sum','shp','fef','ene','tnr','tro','res','ind','slv','agl','agr','ags','awb','swd']
CAMS_sectors_full = ['agl', 'ags', 'awb', 'ene', 'fef', 'ind', 'res', 'shp', 'slv', 'sum', 'swd', 'tnr', 'tro']
#CAMS_sectors_out = ['shp'] #unmask for shp since we will use CAMS shp emissions everywhere
#CAMS_sectors_out = ['tnr','tro']
CAMS_sectors_out = ['agl', 'ags', 'awb', 'ene', 'fef', 'ind', 'res', 'shp', 'slv', 'swd']

CAMS_sectors_names = {'agl':'agriculture_livestock', 'ags':'agriculture_soils', 
                      'awb':'agriculture_wasteburning', 'ene':'power_gen', 
                      'fef':'fugitives', 'ind':'industry', 
                      'res':'residential_commerc_comb', 'shp':'ships', 
                      'slv':'solvents', 'sum':'total', 'swd':'solidwaste_wastewater',
                      'tnr':'offroad', 'tro':'onroad'}

#Settings for which sectors to mask and which not to mask
sector_masking = {'agl':True, 'ags':True, 
                      'awb':True, 'ene':True, 
                      'fef':True, 'ind':True, 
                      'res':True, 'shp':False, 
                      'slv':True, 'sum':True, 'swd':True,
                      'tnr':True, 'tro':True}

# Adds additional suffix to filename for differentiation 
sector_fn_suffix = {'agl':'_masked', 'ags':'_masked', 
                      'awb':'_masked', 'ene':'_masked', 
                      'fef':'_masked', 'ind':'_masked', 
                      'res':'_masked', 'shp':'_unmasked', 
                      'slv':'_masked', 'sum':'_masked', 'swd':'_masked',
                      'tnr':'_masked', 'tro':'_masked'}

# Which sectors should have the mask disabled for CO2 specifically
sector_co2_masking_off = []
sector_ch4_masking_off = ['agl', 'ags']


sector_enable_DOW   = ['tnr','tro']
DOW                = 'sundy'
DOW_diurnal_suffix  = {'sundy': '_FIVE_sundy',  # 
                        'satdy':'_FIVE_satdy',
                        'weekdy':'_FIVE_weekdy'}

# 2020 sector scalings
covid_scalings_fn           = input_dir + 'CONFORM_Glb_0.1x0.1_anthro_adjustment-factors_v1.1_monthly.nc'
covid_scaling_sector_match  = {'ene':'ene_AF_avg', 'ind':'ind_AF_avg',
                               'res':'res_AF_avg','shp':'shp_AF_avg',
                               'tro':'tro_AF_avg'} 
covid_scaling_month_index   = 0

#diurnal_fn         ='CAMS_diurnalweights.csv'
diurnal_fn_base         ='CAMS_diurnalweights' # base of fn, DOW_diurnal_suffix + '.csv' added if enabled for sector
#diurnal_fn        ='CAMS_diurnalweights_FIVE_weekdy.csv'
#diurnal_fn        ='CAMS_diurnalweights_FIVE_satdy.csv'
#diurnal_fn        ='CAMS_diurnalweights_FIVE_sundy.csv'
tz_names = ['Hawaii','Alaska','Pacific','Mountain','Central','Eastern']
tz_offset = {'Hawaii':-10,'Alaska':-9,'Pacific':-8,
             'Mountain':-7,'Central':-6,'Eastern':-5}
tz_mask_fn = 'timezone_mask_nei04k.nc4'

# CONUS Mask File
mask_fn = '/wrk/charkins/emissions/CAMSv4.2/CONUS_mask_NEI04k.nc4'
# Species Conversion Table
species_table_fn = input_dir + 'cams_specieslist_fixedOM.csv'




# Options
mask_on            = True
avg_cams           = False
add_diurnal        = True
new_weights        = False
compress_output    = True
sector_COVIDscalings = False

# If avg_cams = True, specify what months to average over
CAMS_time_start   = -24  # Position of start time of interest (Jan 2019)
CAMS_time_end     = -13   # Position of end time of interest  (Dec 2020)

# If avg_cams = False, specify the position of month of interest
CAMS_time   = -24  # Position of start time of interest (Dec 2020)

mw_for_voc_total = {'HC01':16.,
        'HC02':30.,'HC05':86.,'HC06':44.,'HC07':28.,'HC08':40.,'HC09':56.,'HC10':68.,'HC11':136.,'HC14':30.,'HC18':72.,
        'HC20':184.,'HC31':59.,'HC37':26.,'HC38':78.,'HC39':58.,'HC40':70.,'HC41':92.,'HC42':106.,'HC44':106.,'HC45':106.,
        'HC48':32.,'HC49':81.,'HC50':68.}

# %% Read files
in_ds = xr.open_dataset(CAMS_dir + CAMS_filename_prefix + 'nox' +CAMS_filename_suffix, chunks={'time': 1})
CAMS_nc_struct = in_ds.isel(time=[249]).squeeze(dim='time',drop=True).drop(CAMS_sectors_full)# Make structure for rearranging all of the CAMS info
in_ds.close()

out_example = xr.open_dataset(out_example_fn, chunks={'Time': 1})
times00z = out_example['Times'].values
times12z = xr.open_dataset(out_example_fn2,chunks={'Time':1})['Times'].values
species_table  = pd.read_csv(species_table_fn)  
# grid centers and grid corners for CAMS Domain
in_centers_ds = xr.open_dataset(input_dir + 'gridcenters.nc4')
in_corners_ds = xr.open_dataset(input_dir + 'gridcorners.nc4')

 # %% Make projections between domains

 
#Calculating Projection from WRF LCC to WRF lat long space 
wrf_proj = pyproj.Proj(proj='lcc', # projection type: Lambert Conformal Conic
                       lat_1=out_example.TRUELAT1, lat_2=out_example.TRUELAT2, # Cone intersects with the sphere
                       lat_0=out_example.MOAD_CEN_LAT, lon_0=out_example.STAND_LON, # Center point
                       a=6370000, b=6370000) # Radius
# More info here: https://fabienmaussion.info/2018/01/06/wrf-projection/

cams_proj = pyproj.Proj(proj='latlong',ellps='WGS84',datum ='WGS84') # 

# %% Construct Grid for WRF
#e_spherical, n_spherical = pyproj.transform(cams_proj, wrf_spherical, out_example.CEN_LON, out_example.CEN_LAT)
e, n = pyproj.transform(cams_proj,wrf_proj,out_example.CEN_LON, out_example.CEN_LAT)

# Grid parameters
dx_wrf, dy_wrf = out_example.DX, out_example.DY
nx_wrf, ny_wrf = out_example.dims['west_east'], out_example.dims['south_north']
# Down left corner of the domain
x0_wrf = -(nx_wrf-1) / 2. * dx_wrf + e
y0_wrf = -(ny_wrf-1) / 2. * dy_wrf + n

# Grid of Grid Centers
xx_wrf, yy_wrf = np.meshgrid(np.arange(nx_wrf) * dx_wrf + x0_wrf, np.arange(ny_wrf) * dy_wrf + y0_wrf)
#Transformation of Center X-Y to Center Lat-Lon
#lon_wrf, lat_wrf = wrf_proj(xx_wrf, yy_wrf, inverse=True)
lon_wrf, lat_wrf = pyproj.transform(wrf_proj,cams_proj,xx_wrf,yy_wrf)

# Calculating the boundary X-Y Coordinates
x_b_wrf, y_b_wrf = np.meshgrid(np.arange(nx_wrf+1) * dx_wrf + x0_wrf -dx_wrf/2, np.arange(ny_wrf+1) * dy_wrf + y0_wrf -dy_wrf/2)
#Transformation of Boundary X-Y to Boundary Lat_Lon
lon_b_wrf, lat_b_wrf = pyproj.transform(wrf_proj,cams_proj,x_b_wrf,y_b_wrf)

# %% Construct Grid for CAMS
lon_cams, lat_cams = in_centers_ds['lon_c'].values, in_centers_ds['lat_c'].values
lon_b_cams, lat_b_cams = in_corners_ds['lon_b'].values, in_corners_ds['lat_b'].values

# %% Putting information for regridding in the correct format for xESMF

# Grid Spacing For Output Grid - input for regridder
grid_out = {'lat': lat_wrf, #Center Point Spacing Lat
                    'lon': lon_wrf, #Center Point Spacing Lon
                    'lat_b': lat_b_wrf, # Boundary Spacing Lat 
                    'lon_b': lon_b_wrf, # Boundary Spacing Lon
                    }

# Grid Spacing For Input Grid - input for regridder
grid_in = {'lat': lat_cams, #Center Point Spacing Lat
                    'lon': lon_cams, #Center Point Spacing Lon
                    'lat_b': lat_b_cams, # Boundary Spacing Lat 
                    'lon_b': lon_b_cams, # Boundary Spacing Lon
                    }

# %% Constructing the regridder
if new_weights:
    regridder = xe.Regridder(grid_in, grid_out, method='conservative', reuse_weights=False)
    print("New regrid weights file calculated.")
else:
    regridder = xe.Regridder(grid_in, grid_out, method='conservative', reuse_weights=True)
    print("New regrid weights not calculated, using previous file.")

# %% Construct dataset by sector, conversion to WRF units and species conversion to WRF species

wrf_names = pd.unique(species_table.BIN)

# Loop over each sector and create a sector file with all species
for sector in CAMS_sectors_out:
    out_filename = out_dir + CAMS_sectors_names[sector] + sector_fn_suffix[sector] +out_filename_suffix
    
    #Loop over each species and create each species array within sector
    sector_ds = CAMS_nc_struct.copy()
    
    for species_wrf in WRF_species:
        subtable = species_table.loc[species_table['BIN'] == species_wrf]
        ind = 0
        print('Writing wrf variabe: ' + species_wrf)
        
        # loop over all cams species corresponding to species_wrf
        for species_cams in subtable.cams_name:
            print('Adding emissions from CAMS variable: ' + species_cams)
            
            # Add unit conversion and weighting.
            # Unit convertsion does not include per m^2 to per grid cell, this is done later
            if ind == 0:
                
                unitconv = subtable.loc[subtable['cams_name'] == species_cams].unit_convert_factor.values
                weight = subtable.loc[subtable['cams_name'] == species_cams].weight.values
                print('Weighting: ' + str(weight), '   Unit Conversion Factor: '+ str(unitconv))

                if avg_cams:
                    sector_addition = xr.open_dataset(CAMS_dir + CAMS_filename_prefix + species_cams +CAMS_filename_suffix,
                                                     chunks={'time': 1})[sector][CAMS_time_start:CAMS_time_end,:,:].mean(dim='time',skipna=True)
                elif not avg_cams:
                    sector_addition = xr.open_dataset(CAMS_dir + CAMS_filename_prefix + species_cams +CAMS_filename_suffix,
                                     chunks={'time': 1})[sector][CAMS_time,:,:]#.squeeze(dim='time')
                sector_addition.values = sector_addition.values * weight * unitconv 
                sector_ds[species_wrf] = sector_addition
                sector_addition.close()

                if subtable.cams_name.size > 1:
                    ind = ind + 1
                    continue
            elif ind > 1:
                unitconv = subtable.loc[subtable['cams_name'] == species_cams].unit_convert_factor.values
                weight = subtable.loc[subtable['cams_name'] == species_cams].weight.values
                print('Weighting: ' + str(weight), '   Unit Conversion Factor: '+ str(unitconv))
                if avg_cams:
                    sector_addition = xr.open_dataset(CAMS_dir + CAMS_filename_prefix + species_cams +CAMS_filename_suffix,
                                                     chunks={'time': 1})[sector][CAMS_time_start:CAMS_time_end,:,:].mean(dim='time',skipna=True)
                elif not avg_cams:
                    sector_addition = xr.open_dataset(CAMS_dir + CAMS_filename_prefix + species_cams +CAMS_filename_suffix,
                                     chunks={'time': 1})[sector][CAMS_time,:,:]#.squeeze(dim='time')
                sector_addition.values = sector_addition.values * weight * unitconv
                sector_ds[species_wrf] = sector_ds[species_wrf] + sector_addition
                sector_addition.close()
            #End if
        # End loop over Cams Species
    # End loop over WRF Species
    
    # If 2020 sector scalings are turned on then multiply
    if sector in covid_scaling_sector_match and sector_COVIDscalings:
        COVID_scalings = xr.open_dataset(covid_scalings_fn)
        for species_wrf in WRF_species:
            sector_ds[species_wrf] = (  COVID_scalings[covid_scaling_sector_match[sector]][covid_scaling_month_index,:,:].drop('time') *
                                      sector_ds[species_wrf] )
        COVID_scalings.close()
    # End if

    # Call Regridder on the sector dataset and convert to per grid cell units
    # Assumes grid size of 4km * 4km. This needs to be reversed in the same way
    # when going back to mol km-2 hr-1
    area = 4**2*10**6 # conversion from per m2 to per grid cell
    ds_out = regridder(sector_ds, keep_attrs=True)*area
    sector_ds = sector_ds.drop(WRF_species) #Reset species in sector_ds to reduce memory usage
   
    # If Masking on then mask off everything inside CONUS
    if mask_on:
        #Check if this specific sector should be masked
        if sector_masking[sector]:

            conus_mask = xr.open_dataset(mask_fn)['CONUS_mask'].drop(['lon','lat'])

            # If this sector should not have it's CO2 or CH4 masked then
            if sector in sector_co2_masking_off and sector in sector_ch4_masking_off:
                WRF_species_nomask_removed = list(WRF_species)
                WRF_species_nomask_removed.remove('CO2')
                WRF_species_nomask_removed.remove('HC01')
                
                for species_tomask in WRF_species_nomask_removed :
                    ds_out[species_tomask] = xr.where(conus_mask, 0, ds_out[species_tomask])
            # If this sector should not have it's just CO2 masked then
            elif sector in sector_co2_masking_off and sector not in sector_ch4_masking_off:
                WRF_species_nomask_removed = list(WRF_species)
                WRF_species_nomask_removed.remove('CO2')
                
                for species_tomask in WRF_species_nomask_removed :
                    ds_out[species_tomask] = xr.where(conus_mask, 0, ds_out[species_tomask])
                    
            # If this sector should not have it's just CH4 masked then      
            elif sector not in sector_co2_masking_off and sector in sector_ch4_masking_off:
                WRF_species_nomask_removed = list(WRF_species)
                WRF_species_nomask_removed.remove('HC01')
                
                for species_tomask in WRF_species_nomask_removed:
                    ds_out[species_tomask] = xr.where(conus_mask, 0, ds_out[species_tomask])                    
            else:
                ds_out = xr.where(conus_mask,0,ds_out)
            #end no masking condition if
        #end sector specific masking if
        conus_mask.close()
    #end masking if

    # %% Renaming after regridder to match wrfchemi dimentions and coordinate names
    ds_out = ds_out.rename({'y':'south_north'})
    ds_out = ds_out.rename({'x':'west_east'})
    ds_out = ds_out.rename({'lat':'XLAT'})
    ds_out = ds_out.rename({'lon':'XLONG'})
    ds_out = ds_out.reset_coords(names=('XLAT','XLONG'))
    ds_out = ds_out.transpose('south_north','west_east')
    
        # Move data into 2 12 hour days
    ds_out00z = xr.combine_by_coords([ ds_out.assign_coords({'Time':[0]}), ds_out.assign_coords({'Time':[1]}), 
                                   ds_out.assign_coords({'Time':[2]}), ds_out.assign_coords({'Time':[3]}), 
                                   ds_out.assign_coords({'Time':[4]}),ds_out.assign_coords({'Time':[5]}), 
                                   ds_out.assign_coords({'Time':[6]}),ds_out.assign_coords({'Time':[7]}), 
                                   ds_out.assign_coords({'Time':[8]}),ds_out.assign_coords({'Time':[9]}), 
                                   ds_out.assign_coords({'Time':[10]}), ds_out.assign_coords({'Time':[11]})
                                   ])
    ds_out12z = xr.combine_by_coords([  ds_out.assign_coords({'Time':[12]}), ds_out.assign_coords({'Time':[13]}),
                                   ds_out.assign_coords({'Time':[14]}), ds_out.assign_coords({'Time':[15]}),
                                   ds_out.assign_coords({'Time':[16]}), ds_out.assign_coords({'Time':[17]}),
                                   ds_out.assign_coords({'Time':[18]}), ds_out.assign_coords({'Time':[19]}),
                                   ds_out.assign_coords({'Time':[20]}), ds_out.assign_coords({'Time':[21]}),
                                   ds_out.assign_coords({'Time':[22]}), ds_out.assign_coords({'Time':[23]})
                                   ])
    ds_out = ds_out.drop(WRF_species) #Reset species in ds_out to reduce memory usage
   
    # If diurnal profiling is on, add diurnal profile for sector across each timezone
    if add_diurnal:
        if sector in sector_enable_DOW:
            if DOW == 'weekdy':
                diurnal = pd.read_csv(input_dir + diurnal_fn_base + DOW_diurnal_suffix[DOW] + '.csv') # Import Diurnal Weight File
                print('Diurnal file used: ' + input_dir + diurnal_fn_base + DOW_diurnal_suffix[DOW] + '.csv')
                diurnal = diurnal.append(diurnal)
            elif DOW == 'satdy':
                diurnal_daybefore = pd.read_csv(input_dir + diurnal_fn_base + DOW_diurnal_suffix['weekdy'] + '.csv') # Import Diurnal Weight File
                diurnal = pd.read_csv(input_dir + diurnal_fn_base + DOW_diurnal_suffix[DOW] + '.csv') # Import Diurnal Weight File
                print('Diurnal file 1 used: ' + input_dir + diurnal_fn_base + DOW_diurnal_suffix['weekdy'] + '.csv')
                print('Diurnal file 2 used: ' + input_dir + diurnal_fn_base + DOW_diurnal_suffix[DOW] + '.csv')
                diurnal = diurnal.append(diurnal_daybefore)
            elif DOW == 'sundy':
                diurnal_daybefore = pd.read_csv(input_dir + diurnal_fn_base + DOW_diurnal_suffix['satdy'] + '.csv') # Import Diurnal Weight File
                diurnal = pd.read_csv(input_dir + diurnal_fn_base + DOW_diurnal_suffix[DOW] + '.csv') # Import Diurnal Weight File
                print('Diurnal file 1 used: ' + input_dir + diurnal_fn_base + DOW_diurnal_suffix['satdy'] + '.csv')
                print('Diurnal file 2 used: ' + input_dir + diurnal_fn_base + DOW_diurnal_suffix[DOW] + '.csv')
                diurnal = diurnal.append(diurnal_daybefore)

        else:
            diurnal = pd.read_csv(input_dir + diurnal_fn_base + '.csv') # Import Diurnal Weight File
            print('Diurnal file used: ' + input_dir + diurnal_fn_base + '.csv')
            diurnal = diurnal.append(diurnal)
        # %% Calculate Diurnal Weightings for sector
        diurnal_weights = diurnal[CAMS_sectors_names[sector]].where(diurnal['DOW']=='Average').values
        diurnal_times_local = diurnal['Hour'].where(diurnal['DOW']=='Average').values
    
        #Loop over timezones
        for jj in tz_names:
            start_pos = 23+ tz_offset[jj]
            end_pos = start_pos + 24
            # Move diurnal profile from local time to UTC time
            #diurnal_times_utc = (diurnal_times_local - tz_offset[jj]) % 24
            tz_mask = xr.open_dataset(input_dir + tz_mask_fn)[jj].rename({'y':'south_north','x':'west_east'})
            # Reorder profile times and profiles
            #diurnal_weights_sort = [x for _,x in sorted(zip(diurnal_times_utc,diurnal_weights))]
            diurnal_weights_sort = diurnal_weights[start_pos:end_pos]
            profile = xr.DataArray(diurnal_weights_sort,dims=['Time'],
                   coords = dict( Time=np.arange(0,24) ) )
            
            ds_out00z = xr.where( tz_mask,ds_out00z.transpose('Time','south_north','west_east')*profile[0:12]*24,
                    ds_out00z.transpose('Time','south_north','west_east') )
            ds_out12z = xr.where( tz_mask, ds_out12z.transpose('Time','south_north','west_east')*profile[12:24]*24,
                    ds_out12z.transpose('Time','south_north','west_east') )
            tz_mask.close()
        # End tz loop
        # drop the variables that were added by tz_mask and drop XLAT and XLONG to reset them 
        ds_out00z = ds_out00z.drop(['Time','lat','lon','lat_c','lon_c','XLAT','XLONG'])
        ds_out12z = ds_out12z.drop(['Time','lat','lon','lat_c','lon_c','XLAT','XLONG'])
    else:
        ds_out00z = ds_out00z.drop('Time','XLAT','XLONG')
        ds_out12z = ds_out12z.drop('Time','XLAT','XLONG')

    # Add XLAT and XLONG straight from the template
    ds_out00z['XLAT'] = xr.DataArray(data=out_example['XLAT'].values,dims=['south_north','west_east'])
    ds_out00z['XLONG'] = xr.DataArray(data=out_example['XLONG'].values,dims=['south_north','west_east'])
    ds_out12z['XLAT'] = xr.DataArray(data=out_example['XLAT'].values,dims=['south_north','west_east'])
    ds_out12z['XLONG'] = xr.DataArray(data=out_example['XLONG'].values,dims=['south_north','west_east'])
    ds_out = ds_out.close()  
    # Make sure that dimentions are in correct order 
    ds_out00z = ds_out00z.transpose('Time','south_north','west_east')
    ds_out12z = ds_out12z.transpose('Time','south_north','west_east')
    
    # Set values for times
    ds_out00z['Times'] = xr.DataArray(data=times00z, dims=['Time'])
    ds_out12z['Times'] = xr.DataArray(data=times12z, dims=['Time'])

    # %% Set Storage Settings for output
    chunk_time, chunk_lev, chunk_y, chunk_x =  1, 1, out_example.sizes['south_north'], out_example.sizes['west_east']
    out_coords = 'XLONG XLAT'
    for ii in ds_out00z.data_vars:
        if compress_output:
            ds_out00z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                             'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates': out_coords}
        else:
            ds_out00z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                             '_FillValue': None , 'coordinates': out_coords}
    
    # calculate total VOC and total PM species
    ds_out00z['VOC']=xr.zeros_like(ds_out00z['NOX'])
    for v in mw_for_voc_total.keys():
        ds_out00z['VOC'] =+ ds_out00z[v]*mw_for_voc_total[v]/(10**6) # convert from moles/hour to mt/hour
        
    ds_out00z['PM25-PRI']=xr.zeros_like(ds_out00z['PM04'])
    ds_out00z['PM10-PRI']=xr.zeros_like(ds_out00z['PM04'])
    ds_out00z['PM25-PRI'].values = ds_out00z['PM04'].values + ds_out00z['PM05'].values + ds_out00z['PM06'].values # add together for totals
    ds_out00z['PM10-PRI'].values = ds_out00z['PM04'].values + ds_out00z['PM05'].values + ds_out00z['PM06'].values # add together for totals
    

    ds_out00z['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
    ds_out00z['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
    ds_out00z['Times'].encoding={'char_dim_name':'DateStrLen'}  
    # Set output file global attributes from output example global attributes
    ds_out00z.attrs = out_example.attrs
    # Set variable attributes to match
    for ii in ds_out00z.data_vars:
        ds_out00z[ii].attrs = out_example[ii].attrs
    ds_out00z.attrs['TITLE'] = " Emission Inventory Regridded from CAMSv4.2 0.1x0.1 to WRF lambert conformal domain" ;
    
    # Then do the same for second half of day
    for ii in ds_out12z.data_vars:
        if compress_output:
            ds_out12z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                             'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates': out_coords}
        else:
            ds_out12z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                            '_FillValue': None , 'coordinates': out_coords}

    # calculate total VOC and total PM species
    ds_out12z['VOC']=xr.zeros_like(ds_out12z['NOX'])
    for v in mw_for_voc_total.keys():
        ds_out12z['VOC'] =+ ds_out12z[v]*mw_for_voc_total[v]/(10**6) # convert from moles/hour to mt/hour
        
    ds_out12z['PM25-PRI']=xr.zeros_like(ds_out12z['PM04'])
    ds_out12z['PM10-PRI']=xr.zeros_like(ds_out12z['PM04'])
    ds_out12z['PM25-PRI'].values = ds_out12z['PM04'].values + ds_out12z['PM05'].values + ds_out12z['PM06'].values # add together for totals
    ds_out12z['PM10-PRI'].values = ds_out12z['PM04'].values + ds_out12z['PM05'].values + ds_out12z['PM06'].values # add together for totals
    
    ds_out12z['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
    ds_out12z['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
    ds_out12z['Times'].encoding={'char_dim_name':'DateStrLen'}   
    # Set output file global attributes from output example global attributes
    ds_out12z.attrs = out_example.attrs
    # Set variable attributes to match
    for ii in ds_out12z.data_vars:
        ds_out12z[ii].attrs = out_example[ii].attrs
    ds_out12z.attrs['TITLE'] = " Emission Inventory Regridded from CAMSv4.2 0.1x0.1 to WRF lambert conformal domain"
    
    # %% Save out the Files to NetCDF
    if sector in sector_enable_DOW:
        DOW_fn = '_' + DOW
    else:
        DOW_fn = ''
    ds_out00z.to_netcdf(out_filename+DOW_fn+'_00z'+out_filename_type,format='netCDF4',engine='netcdf4')
    ds_out00z = ds_out00z.drop(WRF_species)
    ds_out12z.to_netcdf(out_filename+DOW_fn+'_12z'+out_filename_type,format='netCDF4',engine='netcdf4')
    ds_out12z = ds_out12z.drop(WRF_species)
    
    print("Finished writing emissions for", sector, "sector", sep=" ")
    
    # End loop over sectors
print('Program Finished')
# End of Program
