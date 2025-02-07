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

class inputs():
    # WRF file for WRF Domain
    #out_example           = xr.open_dataset('wrfchemi_00z_d01.nc', chunks={'Time': 1,'emissions_zdim': 1})
    input_dir             ='/wrk/charkins/emissions/CAMSv6.2/inputs/' 
    out_example_fn        ='/wrk/d2/charkins/templates/with_additional_VCPs_ffCO2/template_new_00to12Z.nc'
    out_example_fn2       ='/wrk/d2/charkins/templates/with_additional_VCPs_ffCO2/template_new_12to24Z.nc'
    #out_dir               = '/wrk/d2/CAMSv4.2/regrid/'
    out_dir               = '/wrk/users/charkins/emissions/CAMSv6.2/regrid_nomask/[YYYY]/Month[MM]/'
    out_filename_suffix   = ''
    out_filename_type     = '.nc4'

    os.chdir(input_dir)
    # CAMS file for data
    #CAMS_dir               = '/wrk/d2/CAMSv4.2/'
    CAMS_dir               = '/wrk/users/charkins/emissions/CAMSv6.2/'
    CAMS_filename_prefix   = 'CAMS-GLOB-ANT_Glb_0.1x0.1_anthro_'
    CAMS_filename_suffix   ='_v6.2_monthly.nc'

    # Species List from CAMS
    CAMS_species = ['bc', 'ch4', 'co2_excl_short-cycle_org_C', 'co2_excl_short-cycle_org_C', 'co', 'nh3', 'nmvoc', 'nox', 'oc', 'so2',
                     'alcohols', 'ethane', 'propane', 'butanes', 'pentanes', 'hexanes', 'ethene', 'propene', 'acetylene', 'isoprene',
                     'monoterpenes', 'other-alkenes-and-alkynes', 'benzene', 'toluene', 'xylene', 'trimethylbenzene', 'other-aromatics', 'esters', 'ethers', 'chlorinated-hydrocarbons',
                     'formaldehyde', 'other-aldehydes', 'total-ketones', 'total-acids', 'other-vocs']
    #CAMS_species = ['nox','ch4','co']
    #WRF_species = ['NOX','CO','CO2']
    WRF_species = ['NOX','CO','ffCO2','CO2','SO2','NH3','PM04','PM05','PM06','HC01',
            'HC02','HC05','HC06','HC07','HC08','HC09','HC10','HC11','HC14','HC18',
            'HC20','HC31','HC37','HC38','HC39','HC40','HC41','HC42','HC44','HC45',
                    'HC48','HC49','HC50']

    # Offroad and onroad (tnr, tro) need to be done separately due to different DOW profiles
    CAMS_sectors_out = ['tnr','tro']
    #CAMS_sectors_out = ['agl', 'ags', 'awb','com', 'ene', 'fef', 'ind', 'ref','res', 'shp', 'slv', 'swd','swl']

    CAMS_sectors_names = {'agl':'agriculture_livestock', 'ags':'agriculture_soils', 
                          'awb':'agriculture_wasteburning','com':'commerc_comb', 'ene':'power_gen', 
                          'fef':'fugitives', 'ind':'industry', 
                          'ref':'refineries','res':'residential_comb', 'shp':'ships', 
                          'slv':'solvents', 'sum':'total', 'swd':'solidwaste_wastewater','swl':'solidwaste_landfill',
                          'tnr':'offroad', 'tro':'onroad'}

    #Settings for which sectors to mask and which not to mask
    sector_masking = {'agl':False, 'ags':False, 
                          'awb':False,'com':False, 'ene':False, 
                          'fef':False, 'ind':False, 'ref':False,
                          'res':False, 'shp':False, 
                          'slv':False, 'sum':False, 'swd':False,'swl':False,
                          'tnr':False, 'tro':False}

    # Adds additional suffix to filename for differentiation 
    sector_fn_suffix = {'agl':'_unmasked', 'ags':'_unmasked', 
                          'awb':'_unmasked','com':'_unmasked', 'ene':'_unmasked', 
                          'fef':'_unmasked', 'ind':'_unmasked', 'ref':'_unmasked',
                          'res':'_unmasked', 'shp':'_unmasked', 
                          'slv':'_unmasked', 'sum':'_unmasked', 'swd':'_unmasked','swl':'_unmasked',
                          'tnr':'_unmasked', 'tro':'_unmasked'}

    # Which sectors should have the mask disabled for CO2 specifically
    sector_co2_masking_off = []
    sector_ch4_masking_off = []


    sector_enable_DOW   = ['tnr','tro']
    DOW                = 'sundy'
    DOW_diurnal_suffix  = {'sundy': '_FIVE_sundy',  # 
                            'satdy':'_FIVE_satdy',
                            'weekdy':'_FIVE_weekdy'}

    # 2020 sector scalings
    #covid_scalings_fn           = input_dir + 'CONFORM_Glb_0.1x0.1_anthro_adjustment-factors_v1.1_monthly.nc'
    #covid_scaling_sector_match  = {'ene':'ene_AF_avg', 'ind':'ind_AF_avg',
    #                               'res':'res_AF_avg','shp':'shp_AF_avg',
    #                               'tro':'tro_AF_avg'} 
    #covid_scaling_month_index   = 0

    #diurnal_fn         ='CAMS_diurnalweights.csv'
    diurnal_fn_base         ='CAMS_diurnalweights' # base of fn, DOW_diurnal_suffix + '.csv' added if enabled for sector
    #diurnal_fn        ='CAMS_diurnalweights_FIVE_weekdy.csv'
    #diurnal_fn        ='CAMS_diurnalweights_FIVE_satdy.csv'
    #diurnal_fn        ='CAMS_diurnalweights_FIVE_sundy.csv'
    tz_names = ['Hawaii','Alaska','Pacific','Mountain','Central','Eastern']
    tz_offset = {'Hawaii':-10,'Alaska':-9,'Pacific':-8,
                 'Mountain':-7,'Central':-6,'Eastern':-5} # daylight savings taken care of in code
    tz_mask_fn = 'timezone_mask_nei04k.nc4'

    # CONUS Mask File
    mask_fn = '/wrk/charkins/emissions/CAMSv4.2/CONUS_mask_NEI04k.nc4'
    # Species Conversion Table
    species_table_fn = input_dir + 'cams_specieslist_v6.2.csv'




    # Options
    mask_on            = False
    #avg_cams           = False
    add_diurnal        = True
    new_weights        = False
    compress_output    = True
    sector_COVIDscalings = False

    # If avg_cams = True, specify what months to average over
    CAMS_time_start   = -24  # Position of start time of interest (Jan 2019)
    CAMS_time_end     = -13   # Position of end time of interest  (Dec 2020)
    
    years = [2020,2021,2022,2023]#[2005,2010,2015,2019]
    months = [1,2,3,4,5,6,7,8,9,10,11,12]
    
    #CAMS_times = {'2005_01':-189,'2005_02':-188,'2006_01':-176,'2006_02':-175,}

    ## If avg_cams = False, specify the position of month of interest
    #CAMS_time   = -189  # Position of start time of interest (Dec 2020)

def construct_regrid():
    
    # %% Read files
    

    out_example = xr.open_dataset(inputs.out_example_fn, chunks={'Time': 1},cache=False)
    
    
    # grid centers and grid corners for CAMS Domain
    in_centers_ds = xr.open_dataset(inputs.input_dir + 'gridcenters.nc4')
    in_corners_ds = xr.open_dataset(inputs.input_dir + 'gridcorners.nc4')

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
    if inputs.new_weights:
        regridder = xe.Regridder(grid_in, grid_out, method='conservative', reuse_weights=False)
        print("New regrid weights file calculated.")
    else:
        regridder = xe.Regridder(grid_in, grid_out, method='conservative', reuse_weights=True)
        print("New regrid weights not calculated, using previous file.")

    return regridder
    
def make_emis(regridder,out_example,times00z,times12z,CAMS_nc_struct,species_table,wrf_names,year,month,sector):
        out_filename = inputs.out_dir.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)) + inputs.CAMS_sectors_names[sector] + inputs.sector_fn_suffix[sector] + inputs.out_filename_suffix
        
        #Loop over each species and create each species array within sector
        sector_ds = CAMS_nc_struct.copy()
        
        for species_wrf in inputs.WRF_species:
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

                    #if inputs.avg_cams:
                    #    sector_addition = xr.open_dataset(inputs.CAMS_dir + inputs.CAMS_filename_prefix + species_cams +inputs.CAMS_filename_suffix,
                    #                                     chunks={'time': 1},cache=False)[sector][inputs.CAMS_time_start:inputs.CAMS_time_end,:,:].mean(dim='time',skipna=True)
                    #sector_addition = xr.open_dataset(inputs.CAMS_dir + inputs.CAMS_filename_prefix + species_cams +inputs.CAMS_filename_suffix,
                    #                 chunks={'time': 1},cache=False)[sector][inputs.CAMS_times[str(year)+'_'+str(month).zfill(2)],:,:]#.squeeze(dim='time')
                    data = xr.open_dataset(inputs.CAMS_dir + inputs.CAMS_filename_prefix + species_cams +inputs.CAMS_filename_suffix,
                                     chunks={'time': 1},cache=False)#[sector].loc[dict(time=str(year)+"-"+str(month).zfill(2)+"-01")]
                    if sector in list(data.data_vars):
                        sector_addition=data[sector].loc[dict(time=str(year)+"-"+str(month).zfill(2)+"-01")]
                    else:
                        sector_addition=data["sum"].loc[dict(time=str(year)+"-"+str(month).zfill(2)+"-01")]
                        sector_addition.values = sector_addition.values*0
                        
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
                    #if avg_cams:
                    #    sector_addition = xr.open_dataset(inputs.CAMS_dir + inputs.CAMS_filename_prefix + species_cams +inputs.CAMS_filename_suffix,
                    #                                     chunks={'time': 1},cache=False)[sector][CAMS_time_start:CAMS_time_end,:,:].mean(dim='time',skipna=True)
                    #elif not avg_cams:
                    data = xr.open_dataset(inputs.CAMS_dir + inputs.CAMS_filename_prefix + species_cams +inputs.CAMS_filename_suffix,
                                     chunks={'time': 1},cache=False)#[sector].loc[dict(time=str(year)+"-"+str(month).zfill(2)+"-01")]#.squeeze(dim='time')
                    if sector in list(data.data_vars):
                        sector_addition=data[sector].loc[dict(time=str(year)+"-"+str(month).zfill(2)+"-01")]
                    else:
                        sector_addition=data["sum"].loc[dict(time=str(year)+"-"+str(month).zfill(2)+"-01")]
                        sector_addition.values = sector_addition.values*0
                        
                    sector_addition.values = sector_addition.values * weight * unitconv
                    sector_ds[species_wrf] = sector_ds[species_wrf] + sector_addition
                    sector_addition.close()
                #End if
            # End loop over Cams Species
        # End loop over WRF Species
        
        # # If 2020 sector scalings are turned on then multiply
        # if sector in inputs.covid_scaling_sector_match and inputs.sector_COVIDscalings:
            # COVID_scalings = xr.open_dataset(inputs.covid_scalings_fn,cache=False)
            # for species_wrf in WRF_species:
                # sector_ds[species_wrf] = (  COVID_scalings[covid_scaling_sector_match[sector]][covid_scaling_month_index,:,:].drop('time') *
                                          # sector_ds[species_wrf] )
            # COVID_scalings.close()
        # # End if

        # Call Regridder on the sector dataset and convert to per grid cell units
        # Assumes grid size of 4km * 4km. This needs to be reversed in the same way
        # when going back to mol km-2 hr-1
        area = 4**2*10**6 # conversion from per m2 to per grid cell
        ds_out = regridder(sector_ds, keep_attrs=True)*area
        sector_ds = sector_ds.drop(inputs.WRF_species) #Reset species in sector_ds to reduce memory usage
       
        # If Masking on then mask off everything inside CONUS
        if inputs.mask_on:
            #Check if this specific sector should be masked
            if inputs.sector_masking[sector]:

                conus_mask = xr.open_dataset(inputs.mask_fn,cache=False)['CONUS_mask'].drop(['lon','lat'])

                # If this sector should not have it's CO2 or CH4 masked then
                if sector in inputs.sector_co2_masking_off and sector in inputs.sector_ch4_masking_off:
                    WRF_species_nomask_removed = list(WRF_species)
                    WRF_species_nomask_removed.remove('CO2')
                    WRF_species_nomask_removed.remove('HC01')
                    
                    for species_tomask in WRF_species_nomask_removed :
                        ds_out[species_tomask] = xr.where(conus_mask, 0, ds_out[species_tomask])
                # If this sector should not have it's just CO2 masked then
                elif sector in inputs.sector_co2_masking_off and sector not in inputs.sector_ch4_masking_off:
                    WRF_species_nomask_removed = list(WRF_species)
                    WRF_species_nomask_removed.remove('CO2')
                    
                    for species_tomask in WRF_species_nomask_removed :
                        ds_out[species_tomask] = xr.where(conus_mask, 0, ds_out[species_tomask])
                        
                # If this sector should not have it's just CH4 masked then      
                elif sector not in inputs.sector_co2_masking_off and sector in inputs.sector_ch4_masking_off:
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
        ds_out = ds_out.drop(inputs.WRF_species) #Reset species in ds_out to reduce memory usage
       
        # If diurnal profiling is on, add diurnal profile for sector across each timezone
        if inputs.add_diurnal:
            if sector in inputs.sector_enable_DOW:
                if inputs.DOW == 'weekdy':
                    diurnal = pd.read_csv(inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix[inputs.DOW] + '.csv') # Import Diurnal Weight File
                    print('Diurnal file used: ' + inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix[inputs.DOW] + '.csv')
                    diurnal = diurnal.append(diurnal)
                elif inputs.DOW == 'satdy':
                    diurnal_daybefore = pd.read_csv(inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix['weekdy'] + '.csv') # Import Diurnal Weight File
                    diurnal = pd.read_csv(inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix[inputs.DOW] + '.csv') # Import Diurnal Weight File
                    print('Diurnal file 1 used: ' + inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix['weekdy'] + '.csv')
                    print('Diurnal file 2 used: ' + inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix[inputs.DOW] + '.csv')
                    diurnal = diurnal.append(diurnal_daybefore)
                elif inputs.DOW == 'sundy':
                    diurnal_daybefore = pd.read_csv(inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix['satdy'] + '.csv') # Import Diurnal Weight File
                    diurnal = pd.read_csv(inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix[inputs.DOW] + '.csv') # Import Diurnal Weight File
                    print('Diurnal file 1 used: ' + inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix['satdy'] + '.csv')
                    print('Diurnal file 2 used: ' + inputs.input_dir + inputs.diurnal_fn_base + inputs.DOW_diurnal_suffix[inputs.DOW] + '.csv')
                    diurnal = diurnal.append(diurnal_daybefore)

            else:
                diurnal = pd.read_csv(inputs.input_dir + inputs.diurnal_fn_base + '.csv') # Import Diurnal Weight File
                print('Diurnal file used: ' + inputs.input_dir + inputs.diurnal_fn_base + '.csv')
                diurnal = diurnal.append(diurnal)
            # %% Calculate Diurnal Weightings for sector
            diurnal_weights = diurnal[inputs.CAMS_sectors_names[sector]].where(diurnal['DOW']=='Average').values
            diurnal_times_local = diurnal['Hour'].where(diurnal['DOW']=='Average').values
        
            #Loop over timezones
            for jj in inputs.tz_names:
                if 3 <= month <=10:
                    dst_adjust=1
                else:
                    dst_adjust=0
                
                start_pos = 23+ (inputs.tz_offset[jj] +dst_adjust)
                end_pos = start_pos + 24
                # Move diurnal profile from local time to UTC time
                #diurnal_times_utc = (diurnal_times_local - tz_offset[jj]) % 24
                tz_mask = xr.open_dataset(inputs.input_dir + inputs.tz_mask_fn)[jj].rename({'y':'south_north','x':'west_east'})
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
            if inputs.compress_output:
                ds_out00z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                                 'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates': out_coords}
            else:
                ds_out00z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                                 '_FillValue': None , 'coordinates': out_coords}


        ds_out00z['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
        ds_out00z['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
        ds_out00z['Times'].encoding={'char_dim_name':'DateStrLen'}  
        # Set output file global attributes from output example global attributes
        ds_out00z.attrs = out_example.attrs
        # Set variable attributes to match
        for ii in ds_out00z.data_vars:
            ds_out00z[ii].attrs = out_example[ii].attrs
        ds_out00z.attrs['TITLE'] = " Emission Inventory Regridded from CAMSv6.2 0.1x0.1 to WRF lambert conformal domain" ;
        
        # Then do the same for second half of day
        for ii in ds_out12z.data_vars:
            if inputs.compress_output:
                ds_out12z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                                 'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates': out_coords}
            else:
                ds_out12z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                                '_FillValue': None , 'coordinates': out_coords}


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
        if sector in inputs.sector_enable_DOW:
            DOW_fn = '_' + inputs.DOW
        else:
            DOW_fn = ''
        ds_out00z.to_netcdf(out_filename+DOW_fn+'_00z'+inputs.out_filename_type,format='netCDF4',engine='netcdf4')
        ds_out00z = ds_out00z.drop(inputs.WRF_species)
        ds_out12z.to_netcdf(out_filename+DOW_fn+'_12z'+inputs.out_filename_type,format='netCDF4',engine='netcdf4')
        ds_out12z = ds_out12z.drop(inputs.WRF_species)
        
        print("Finished writing emissions for", sector, "sector", sep=" ")
        
    # End loop over sectors
    
def main():
    
    in_ds = xr.open_dataset(inputs.CAMS_dir + inputs.CAMS_filename_prefix + 'nox' +inputs.CAMS_filename_suffix, chunks={'time': 1},cache=False)
    in_ds.close()
    
    
    regridder = construct_regrid()
    
    out_example = xr.open_dataset(inputs.out_example_fn, chunks={'Time': 1},cache=False)
    times00z = out_example['Times'].values
    times12z = xr.open_dataset(inputs.out_example_fn2,chunks={'Time':1})['Times'].values

    # %% Construct dataset by sector, conversion to WRF units and species conversion to WRF species
    CAMS_nc_struct = in_ds.isel(time=[0]).squeeze(dim='time',drop=True).drop(in_ds.data_vars)# Make structure for rearranging all of the CAMS info
    species_table  = pd.read_csv(inputs.species_table_fn)
    wrf_names = pd.unique(species_table.BIN)
    
    for year in inputs.years:
        
        for month in inputs.months:
            
            print('Working on emissions for: '+ str(year) + str(month).zfill(2))
            
            # Loop over each sector and create a sector file with all species
            for sector in inputs.CAMS_sectors_out:
                make_emis(regridder,out_example,times00z,times12z,CAMS_nc_struct,species_table,wrf_names,year,month,sector)
            
    print('Program Finished')
    # End of Program
# End Main

if __name__ == "__main__":
    main()