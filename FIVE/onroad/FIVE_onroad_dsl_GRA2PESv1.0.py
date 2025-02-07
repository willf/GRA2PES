# Uses geopandas_env
import xarray as xr
import numpy as np
import pandas as pd
# import datetime as dt
import time
import gc

class inputs():
    # Input/Output Dirs
    input_dir = '/wrk/charkins/FIVE/FIVE_onroad_python/inputs/'
    output_base_dir = '/wrk/users/charkins/emissions/FIVE/onroad_dsl_ffCO2/'
    # Domain
    domain_dir   = '/wrk/charkins/FIVE/Domain/'
    domain_name  = "nei04k_domain"
    # Template
    template_1   ='/wrk/charkins/FIVE/FIVE_onroad_python/template/with_CO2ff/onroad_00to12Z.nc'
    template_2   ='/wrk/charkins/FIVE/FIVE_onroad_python/template/with_CO2ff/onroad_12to24Z.nc'
    
    #Shape
    five_shp_fn  = '2018_fuel'
    # five_fuel_fn = 'five2019.rds'
    
    # Script will convert dbf to csv because it's faster to import. if option below is 
    # false, it will only write a new csv if one doesn't exist in the input folder
    # if true it will replace the old .csv with a new one loaded from the shapefile dbf
    new_five_shp_csv = False
    
    # Annual Fuel Scaling, this controls the year that monthly fuel scalings are based from
    ann_fuelscale_fn_dict    = {2005:'onroad_fuel2005.csv',2015:'onroad_fuel2018.csv',2016:'onroad_fuel2018.csv',2017:'onroad_fuel2018.csv',2018:'onroad_fuel2018.csv',2019:'onroad_fuel2018.csv',
                                2020:'onroad_fuel2018.csv',2021:'onroad_fuel2018.csv',2022:'onroad_fuel2022.csv'}
    
    # Emission Factors
    ld_ef_fn_base = 'onroad_gasefYYYY_MM_ffCO2.csv'
    ld_dsl_ef_fn_base = 'onroad_lddslefYYYY_MM_ffCO2.csv'
    hd_ef_fn_base = 'onroad_dslefYYYY_ffCO2.csv' 
    
    # Speciation files (number of speciation profiles must match number of sectors)
    voc_wt_file    = 'onroad_wt.csv'
    voc_mw_file    = 'onroad_mw.csv'
    pm25_wt_file   = 'onroad_pm25.csv'
    
    # Diurnal Files (number of diurnal profiles must match number of sectors)
    monthly_fn = 'onroad_monthly.csv'
    monthly_fuel_fn_dict             = {2005:'onroad_monthly_fuelscale_2005to2005.csv',2015:'onroad_monthly_fuelscale_2018to2015.csv',2017:'onroad_monthly_fuelscale_2018to2017.csv',2018:'onroad_monthly_fuelscale_2018to2018.csv', 2019:'onroad_monthly_fuelscale_2018to2019.csv',
                                        2020:'onroad_monthly_fuelscale_2018to2020.csv', 2021:'onroad_monthly_fuelscale_2018to2021.csv', 2022:'onroad_monthly_fuelscale_2022to2022.csv'}
    # monthly_PADDtoSTATE_fn_base = 'onroad_PADDtoSTATE_YYYY.csv'
    monthly_PADDtoSTATE_fn_dict = {2005:'onroad_PADDtoSTATE_2005to2005.csv',2015:'onroad_PADDtoSTATE_2018to2015.csv',2016:'onroad_PADDtoSTATE_2018to2016.csv',2017:'onroad_PADDtoSTATE_2018to2017.csv',2018:'onroad_PADDtoSTATE_2018to2018.csv',2019:'onroad_PADDtoSTATE_2018to2019.csv',
                                   2020:'onroad_PADDtoSTATE_2018to2020.csv', 2021:'onroad_PADDtoSTATE_2018to2021.csv', 2022:'onroad_PADDtoSTATE_2022to2022.csv'}
    diurnal_fn             = 'onroad_diurnal.csv' #local time
    diurnal_march_fn             = 'onroad_diurnal_plus_halfhour_offset_formarch.csv' #local time
    
    years = [2021]
    months = [1,2,3,4,5,6,7,8,9,10,11,12]
    
    # Takes the year and month and gives back file names for emission factor and scale factor file names
    def get_input_filenames(self,year,month):
        hd_ef_fn = self.hd_ef_fn_base.replace('YYYY',str(year))
        ld_ef_fn = self.ld_ef_fn_base.replace('YYYY',str(year)).replace('MM',str(month).zfill(2))
        if self.ld_diesel_on:
            ld_dsl_ef_fn = self.ld_dsl_ef_fn_base.replace('YYYY',str(year)).replace('MM',str(month).zfill(2))
        monthly_fuel_fn  = self.monthly_fuel_fn_dict[year]
        monthly_PADDtoSTATE_fn = self.monthly_PADDtoSTATE_fn_dict[year]
        ann_fuelscale_fn = self.ann_fuelscale_fn_dict[year]
        if self.ld_diesel_on:
            ef_sf_fn_dict = {'hd_ef_fn':hd_ef_fn,'ld_ef_fn':ld_ef_fn,'monthly_fuel_fn':monthly_fuel_fn,
                             'monthly_PADDtoSTATE_fn':monthly_PADDtoSTATE_fn,'ann_fuelscale_fn':ann_fuelscale_fn,
                             'ld_dsl_ef_fn':ld_dsl_ef_fn}
        else:
            ef_sf_fn_dict = {'hd_ef_fn':hd_ef_fn,'ld_ef_fn':ld_ef_fn,'monthly_fuel_fn':monthly_fuel_fn,
                             'monthly_PADDtoSTATE_fn':monthly_PADDtoSTATE_fn,'ann_fuelscale_fn':ann_fuelscale_fn }
        return ef_sf_fn_dict
    
    # Output directory structure, this must already exist. script wont create new variables
    out_struct = {'year':{'name': dict(zip(years, [str(y) for y in years])),
                          'values': years},
                  'month':{'name': dict(zip(months,[('Month' + str(m).zfill(2)) for m in months])),
                          'values': months}}
    out_fn1 = 'onroad_dsl_00to12Z.nc' # File name for the first half of the day
    out_fn2 = 'onroad_dsl_12to24Z.nc' # File name for the second half of the day
    
    # DOW list
    dow_list   = ['weekdy','satdy','sundy']
    day_map = {'weekdy':0,'satdy':1,'sundy':2} #Index corresponding to each day of week in ds, can leave unchanged
    #Timezone info
    tz_offset = {'Hawaii':-10,'Alaska':-9,'Pacific':-8,
             'Mountain':-7,'Central':-6,'Eastern':-5,'Atlantic':-4} # standard time, daylight savings offset added later
    daylight_savings = {1:False,2:False,3:True,4:True,5:True,6:True,
                        7:True,8:True,9:True,10:True,11:False,12:False}
    
    # Flags
    output_on = True
    ld_diesel_on = False
    output_fuel_types = ['bus','dsl'] # Which fuel types should be summed into the output
    dsl_ef_national = True # If diesel EFs are national then don't loop over all states
    dsl_ef_national_type_list = ['dsl','bus'] # which fuel types are national
    all_species = True # Should all species available be output
     # if not all species then list which, all species will be in netcdf output still but other species will have zero value everywhere
    # species_sub_list = ['NOX','CO','NH3','PM10-PRI','PM25-PRI','SO2','VOC']
    species_sub_list = ['NOX','CO','CO2']
    # species_sub_list = ['NOX','CO','CO2','VOC','HC01','PM25-PRI','PM10-PRI']
    low_mem_on = True # when on, templates opened with xr.open_dataset() rather than the xr.load_dataset()
    
    # Options for saving fuel consumption netcdf files and scaling factor files
    save_monthscaled_fuel = True # save fuel use averaged for the month to its directory
    monthscaled_fuel_fn = 'fuel_avg.nc'
    save_diurnal_fuel = False # save fuel use for the month with diurnal scalings applied
    diurnal_fuel_fn = 'fuel_full.nc'
    save_fuel_scaling = True  # save overall fuel scaling for the month to its directory
    fuel_scaling_fn = 'fuel_scalings.csv'
    
# end inputs class

#

# initialize_fuel_ds
# Place shp fuel data in an xarray dataset and convert to kg/d
# No scaling factors are applied here
def initialize_fuel_ds():
    fuel_ds = xr.open_dataset(inputs.template_1).drop_dims(drop_dims=['Time'])
    # Make arrays for fuel, initialized to zero everywhere
    fuel_ds['gas'] = xr.DataArray(data=np.zeros(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east'])),
        dims=['south_north','west_east'],
        attrs=dict(description='GAS',units='kg/d'))
    fuel_ds['dsl'] = xr.DataArray(data=np.zeros(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east'])),
        dims=['south_north','west_east'],
        attrs=dict(description='DSL',units='kg/d'))
    if inputs.ld_diesel_on:
        fuel_ds['ld_dsl'] = xr.DataArray(data=np.zeros(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east'])),
                                      dims=['south_north','west_east'],
                                      attrs=dict(description='LD_DSL',units='kg/d'))
    
    fuel_ds['bus'] = xr.DataArray(data=np.zeros(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east'])),
        dims=['south_north','west_east'],
        attrs=dict(description='BUS',units='kg/d'))
    # fuel_ds['tz'] = xr.DataArray(data=np.empty(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east']),dtype=object),
    #     dims=['south_north','west_east'],
    #     attrs=dict(description='Timezone'))
    fuel_ds['fips'] = xr.DataArray(data=np.zeros(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east'])),
        dims=['south_north','west_east'],
        attrs=dict(description='State FIPS Code'))
    # fuel_ds['urb_rur'] = xr.DataArray(data=np.empty(shape=(fuel_ds.dims['south_north'],fuel_ds.dims['west_east']),dtype=object),
    #     dims=['south_north','west_east'],
    #     attrs=dict(description='Urban or Rural'))
    
    fuel_df = read_shp_tbl()
    
    # Find indicies of interest from df
    ind_x = (fuel_df['col'].values.astype(int)).tolist()
    ind_y = (fuel_df['row'].values.astype(int)).tolist()
    
    # Assign values from df and convert from gal/d to kg/d
    fuel_ds['gas'].data[ind_y,ind_x] = fuel_df['gas'].astype(float).values * 3.785 * 0.74
    fuel_ds['dsl'].data[ind_y,ind_x] = fuel_df['dsl'].astype(float).values * 3.785 * 0.85
    fuel_ds['bus'].data[ind_y,ind_x] = fuel_df['dsl'].astype(float).values * 0
    if inputs.ld_diesel_on:
        fuel_ds['ld_dsl'].data[ind_y,ind_x] = fuel_df['ld_dsl'].astype(float).values * 3.785 * 0.85
    
    # # Assign values for timezone, fips and urb_rur to ds
    # fuel_ds['tz'].data[ind_y,ind_x] = fuel_df['tz'].astype(str).values
    fuel_ds['fips'].data[ind_y,ind_x] = fuel_df['state_fips'].astype(int).values
    # fuel_ds['urb_rur'].data[ind_y,ind_x] = fuel_df['urb_rur'].astype(str).values
    
    # for ii in ['gas','bus','dsl']:
    #     fuel_ds[ii].encoding={'dtype': 'float32', 'chunksizes':(fuel_ds.sizes['south_north'], fuel_ds.sizes['west_east']),
    #                       'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
    # for ii in ['tz','fips','urb_rur']:
    #     fuel_ds[ii].encoding={'chunksizes':(fuel_ds.sizes['south_north'], fuel_ds.sizes['west_east']),
    #                       'zlib': True, 'complevel': 1, '_FillValue': None}
    # fuel_ds.to_netcdf(inputs.output_base_dir+'fuel_fromshp.nc',format='netCDF4',engine='netcdf4')

    return fuel_ds, fuel_df

# end initialize_fuel_ds

# read_shp_tbl
# Reads the .dbf fuel file and returns as a pandas dataframe
def read_shp_tbl():
    from os.path import exists
    exist = exists(inputs.input_dir+inputs.five_shp_fn+'_dbf.csv')
    
    if exist and not inputs.new_five_shp_csv:
        df = pd.read_csv(inputs.input_dir+inputs.five_shp_fn+'_dbf.csv')
    else:
        import geopandas as gpd
        tbl_fn = inputs.input_dir + inputs.five_shp_fn + '.dbf'
        gdf = gpd.read_file(tbl_fn)
        df = pd.DataFrame(gdf).drop(columns=['geometry']) 
        df['gas'] = df['gas'].fillna(value=0)
        df['dsl'] = df['dsl'].fillna(value=0)
        if inputs.ld_diesel_on:
            df['ld_dsl'] = df['ld_dsl'].fillna(value=0)
        #df['bus'] = df['bus'].fillna(value=0)
        df.to_csv(inputs.input_dir+inputs.five_shp_fn+'_dbf.csv',index=False)
    return df
# End read shape tbl

# make_diurnalfactors
# makes diurnal factorsfor each category
def make_diurnalfactors(month):
    if month ==3:  
        diurnal_scale = pd.read_csv(inputs.input_dir+inputs.diurnal_march_fn,
                                    dtype={'Hour':int})
    else:
        diurnal_scale = pd.read_csv(inputs.input_dir+inputs.diurnal_fn,
                                    dtype={'Hour':int})
    
    df_out = pd.DataFrame({'Timezone' : [],
                           'Day' : [],
                           'Hour_UTC' : [],
                           'GAS_URB':[],
                           'GAS_RUR':[],
                           'DSL_URB':[],
                           'DSL_RUR':[],
                           'BUS_URB':[],
                           'BUS_RUR':[]})
    if inputs.ld_diesel_on:
        df_out = pd.DataFrame({'Timezone' : [],
                           'Day' : [],
                           'Hour_UTC' : [],
                           'GAS_URB':[],
                           'GAS_RUR':[],
                           'DSL_URB':[],
                           'DSL_RUR':[],
                           'LD_DSL_URB':[],
                           'LD_DSL_RUR':[],
                           'BUS_URB':[],
                           'BUS_RUR':[]})
    else:
        df_out = pd.DataFrame({'Timezone' : [],
                           'Day' : [],
                           'Hour_UTC' : [],
                           'GAS_URB':[],
                           'GAS_RUR':[],
                           'DSL_URB':[],
                           'DSL_RUR':[],
                           'BUS_URB':[],
                           'BUS_RUR':[]})
    
    # Loop over timezones and make a table for each
    for tz in inputs.tz_offset.keys():
        # If during daylight savings, add +1 hour to utc offset
        if inputs.daylight_savings[month]:
            dst_offset = 1
        else:
            dst_offset = 0
        offset = inputs.tz_offset[tz] + dst_offset
        for day in inputs.dow_list:
            if day =='weekdy':
                df = diurnal_scale.loc[diurnal_scale['Day']=='weekdy',:].reset_index(drop=True)
                df.loc[:,'Hour_UTC'] = (df.loc[:,'Hour']-offset)%24 # Modulo operator to wrap around
                df = df.sort_values(by='Hour_UTC').reset_index(drop=True)
            elif day =='satdy':
                # Cases depending on value of offset
                if offset < 0:
                    df1 = diurnal_scale.loc[diurnal_scale['Day']=='weekdy',:].reset_index(drop=True)
                    df2 = diurnal_scale.loc[diurnal_scale['Day']=='satdy',:].reset_index(drop=True)
                    df1.loc[:,'Hour_UTC'] = (df1['Hour']-offset)%24 # Modulo operator to wrap around
                    df1 = df1.sort_values(by='Hour_UTC').reset_index(drop=True)
                    df2.loc[:,'Hour_UTC'] = (df2['Hour']-offset)%24 # Modulo operator to wrap around
                    df2 = df2.sort_values(by='Hour_UTC').reset_index(drop=True)
                    # df = pd.concat([df1[0:np.abs(offset),:],df2[np.abs(offset):24,:]])
                    df = pd.concat([df1.head(np.abs(offset)),df2.tail(24-np.abs(offset))])
                elif offset == 0:
                    df = diurnal_scale.loc[diurnal_scale['Day']=='satdy',:].reset_index(drop=True)
                    df.loc[:,'Hour_UTC'] = (df.loc[:,'Hour']-offset)%24 # Modulo operator to wrap around
                    df = df.sort_values(by='Hour_UTC').reset_index(drop=True)
                elif offset > 0:
                    df1 = diurnal_scale.loc[diurnal_scale['Day']=='satdy',:].reset_index(drop=True)
                    df2 = diurnal_scale.loc[diurnal_scale['Day']=='sundy',:].reset_index(drop=True)
                    df1.loc[:,'Hour_UTC'] = (df1.loc[:,'Hour']-offset)%24 # Modulo operator to wrap around
                    df1 = df1.sort_values(by='Hour_UTC').reset_index(drop=True)
                    df2.loc[:,'Hour_UTC'] = (df2.loc[:,'Hour']-offset)%24 # Modulo operator to wrap around
                    df2 = df2.sort_values(by='Hour_UTC').reset_index(drop=True)
                    # df = pd.concat([df1[np.abs(offset):24,:],df2[0:np.abs(offset),:]])
                    df = pd.concat([df1.tail(24-np.abs(offset)),df2.head(np.abs(offset))]).reset_index(drop=True)
                #end if
            elif day =='sundy':
                # Cases depending on value of offset
                if offset < 0:
                    df1 = diurnal_scale.loc[diurnal_scale['Day']=='satdy',:].reset_index(drop=True)
                    df2 = diurnal_scale.loc[diurnal_scale['Day']=='sundy',:].reset_index(drop=True)
                    df1.loc[:,'Hour_UTC'] = (df1['Hour']-offset)%24 # Modulo operator to wrap around
                    df1 = df1.sort_values(by='Hour_UTC').reset_index(drop=True)
                    df2.loc[:,'Hour_UTC'] = (df2['Hour']-offset)%24 # Modulo operator to wrap around
                    df2 = df2.sort_values(by='Hour_UTC').reset_index(drop=True)
                    # df = pd.concat([df1[0:np.abs(offset),:],df2[np.abs(offset):24,:]])
                    df = pd.concat([df1.head(np.abs(offset)),df2.tail(24-np.abs(offset))]).reset_index(drop=True)
                elif offset == 0:
                    df = diurnal_scale.loc[diurnal_scale['Day']=='sundy',:].reset_index(drop=True)
                    df.loc[:,'Hour_UTC'] = (df['Hour']-offset)%24 # Modulo operator to wrap around
                    df = df.sort_values(by='Hour_UTC').reset_index(drop=True)
                elif offset > 0:
                    df1 = diurnal_scale.loc[diurnal_scale['Day']=='sundy',:].reset_index(drop=True)
                    df2 = diurnal_scale.loc[diurnal_scale['Day']=='weekdy',:].reset_index(drop=True)
                    df1.loc[:,'Hour_UTC'] = (df1['Hour']-offset)%24 # Modulo operator to wrap around
                    df1 = df1.sort_values(by='Hour_UTC').reset_index(drop=True)
                    df2.loc[:,'Hour_UTC'] = (df2['Hour']-offset)%24 # Modulo operator to wrap around
                    df2 = df2.sort_values(by='Hour_UTC').reset_index(drop=True)
                    # df = pd.concat([df1[np.abs(offset):24,:],df2[0:np.abs(offset),:]])
                    df = pd.concat([df1.tail(24-np.abs(offset)),df2.head(np.abs(offset))]).reset_index(drop=True)
                #end if
            #end if
            if inputs.ld_diesel_on:
                df_append = pd.DataFrame({'Timezone' : [tz]*24,
                           'Day' : [day]*24,
                           'Hour_UTC' : df['Hour_UTC'].values,
                           'GAS_URB':df['GAS_URB'].values,
                           'GAS_RUR':df['GAS_RUR'].values,
                           'DSL_URB':df['DSL_URB'].values,
                           'DSL_RUR':df['DSL_RUR'].values,
                           'LD_DSL_URB':df['LD_DSL_URB'].values,
                           'LD_DSL_RUR':df['LD_DSL_RUR'].values,
                           'BUS_URB':df['BUS_URB'].values,
                           'BUS_RUR':df['BUS_RUR'].values})
            else:
                df_append = pd.DataFrame({'Timezone' : [tz]*24,
                           'Day' : [day]*24,
                           'Hour_UTC' : df['Hour_UTC'].values,
                           'GAS_URB':df['GAS_URB'].values,
                           'GAS_RUR':df['GAS_RUR'].values,
                           'DSL_URB':df['DSL_URB'].values,
                           'DSL_RUR':df['DSL_RUR'].values,
                           'BUS_URB':df['BUS_URB'].values,
                           'BUS_RUR':df['BUS_RUR'].values})
            df_out = pd.concat([df_out,df_append])
            del df
        # end for over days
    # end for over timezones
    # print(df_out)
    return df_out
# end make_diurnalfactors

# make_emissionfactors
# makes emission factors for each category
def make_emissionfactors(year,month):
    
    # Get the file names needed
    ef_sf_fn_dict = inputs().get_input_filenames(year = year,month = month)
    # Read in efs as dataframes
    ld_df = pd.read_csv(inputs.input_dir+ef_sf_fn_dict['ld_ef_fn'])
    hd_df = pd.read_csv(inputs.input_dir+ef_sf_fn_dict['hd_ef_fn'])
    if inputs.ld_diesel_on:
        ld_dsl_df = pd.read_csv(inputs.input_dir+ef_sf_fn_dict['ld_dsl_ef_fn'])
    
    voc_wt = pd.read_csv(inputs.input_dir+inputs.voc_wt_file).fillna(0)
    voc_mw = pd.read_csv(inputs.input_dir+inputs.voc_mw_file).fillna(0)
    pm25_wt = pd.read_csv(inputs.input_dir+inputs.pm25_wt_file).fillna(0)
    
    # Combine into summed category PM25 and PM10
    ld_df.loc[:,'PM25-PRI'] = ld_df.loc[:,'PM25_PRI_EXH'] + ld_df.loc[:,'PM25_PRI_TIRE'] + ld_df.loc[:,'PM25_PRI_BRAKE']
    ld_df.loc[:,'PM10-PRI'] = ld_df.loc[:,'PM10_PRI_EXH'] + ld_df.loc[:,'PM10_PRI_TIRE'] + ld_df.loc[:,'PM10_PRI_BRAKE']
    hd_df.loc[:,'PM25-PRI'] = hd_df.loc[:,'PM25_PRI_EXH'] + hd_df.loc[:,'PM25_PRI_TIRE'] + hd_df.loc[:,'PM25_PRI_BRAKE']
    hd_df.loc[:,'PM10-PRI'] = hd_df.loc[:,'PM10_PRI_EXH'] + hd_df.loc[:,'PM10_PRI_TIRE'] + hd_df.loc[:,'PM10_PRI_BRAKE']
    if inputs.ld_diesel_on:
        ld_dsl_df.loc[:,'PM25-PRI'] = ld_dsl_df.loc[:,'PM25_PRI_EXH'] + ld_dsl_df.loc[:,'PM25_PRI_TIRE'] + ld_dsl_df.loc[:,'PM25_PRI_BRAKE']
        ld_dsl_df.loc[:,'PM10-PRI'] = ld_dsl_df.loc[:,'PM10_PRI_EXH'] + ld_dsl_df.loc[:,'PM10_PRI_TIRE'] + ld_dsl_df.loc[:,'PM10_PRI_BRAKE']
    
    # Combine into summed category for VOC
    ld_df.loc[:,'VOC'] = ld_df.loc[:,'VOC_EXH'] + ld_df.loc[:,'VOC_EVAP']
    hd_df.loc[:,'VOC'] = hd_df.loc[:,'VOC_EXH'] + hd_df.loc[:,'VOC_EVAP']
    if inputs.ld_diesel_on:
        ld_dsl_df.loc[:,'VOC'] = ld_dsl_df.loc[:,'VOC_EXH'] + ld_dsl_df.loc[:,'VOC_EVAP']
    
    ld_df = ld_df.set_index('FIPS')
    hd_df = hd_df.set_index('FIPS')
    if inputs.ld_diesel_on:
        ld_dsl_df = ld_dsl_df.set_index('FIPS')
    gas_ef_df = ld_df.copy()
    dsl_ef_df = hd_df.copy()
    if inputs.ld_diesel_on:
        ld_dsl_ef_df = ld_dsl_df.copy()
    bus_ef_df = hd_df.copy()
    
    # Loop over PM25 bins to speciate them into g/kg emission factors
    for pmbin in pm25_wt.loc[:,'Bin']:
        gas_ef_df.loc[:,pmbin] = (ld_df.loc[:,'PM25_PRI_EXH']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'GAS_EXH'].values/100 + 
                        ld_df.loc[:,'PM25_PRI_TIRE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'TIRE'].values/100 + 
                        ld_df.loc[:,'PM25_PRI_BRAKE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'BRAKE'].values/100 ) 
        dsl_ef_df.loc[:,pmbin] = (hd_df.loc[:,'PM25_PRI_EXH']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'DSL_EXH'].values/100 + 
                        hd_df.loc[:,'PM25_PRI_TIRE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'TIRE'].values/100 + 
                        hd_df.loc[:,'PM25_PRI_BRAKE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'BRAKE'].values/100 )
        bus_ef_df.loc[:,pmbin] = (hd_df.loc[:,'PM25_PRI_EXH']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'DSL_EXH'].values/100 + 
                        hd_df.loc[:,'PM25_PRI_TIRE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'TIRE'].values/100 + 
                        hd_df.loc[:,'PM25_PRI_BRAKE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'BRAKE'].values/100 )
        if inputs.ld_diesel_on:
            ld_dsl_ef_df.loc[:,pmbin] = (ld_dsl_df.loc[:,'PM25_PRI_EXH']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'LD_DSL_EXH'].values/100 + 
                        ld_dsl_df.loc[:,'PM25_PRI_TIRE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'TIRE'].values/100 + 
                        ld_dsl_df.loc[:,'PM25_PRI_BRAKE']*pm25_wt.loc[pm25_wt['Bin'] == pmbin,'BRAKE'].values/100 ) 
            
        
    # Loop over VOC bins to speciate them into mol/kg emission factors
    for vocbin in voc_wt.loc[:,'Bin']:
        exh = ld_df.loc[:,'VOC_EXH']*voc_wt.loc[voc_wt['Bin'] == vocbin,'GAS_EXH'].values/100/voc_mw.loc[voc_mw['Bin'] == vocbin,'GAS_EXH'].values
        evap = ld_df.loc[:,'VOC_EVAP']*voc_wt.loc[voc_wt['Bin'] == vocbin,'GAS_EVAP'].values/100/voc_mw.loc[voc_mw['Bin'] == vocbin,'GAS_EVAP'].values
        if exh.isna().any():
            exh = exh.fillna(0)
        if evap.isna().any():
            evap=evap.fillna(0)
        gas_ef_df.loc[:,vocbin] = ( evap+exh )
        # VOC_EVAP not included for dsl and bus because there are no weights for it. Evap emission factors for dsl also zero
        dsl_ef_df.loc[:,vocbin] = (hd_df.loc[:,'VOC_EXH']*voc_wt.loc[voc_wt['Bin'] == vocbin,'DSL_EXH'].values/100/voc_mw.loc[voc_mw['Bin'] == vocbin,'DSL_EXH'].values )
        bus_ef_df.loc[:,vocbin] = (hd_df.loc[:,'VOC_EXH']*voc_wt.loc[voc_wt['Bin'] == vocbin,'BUS_EXH'].values/100/voc_mw.loc[voc_mw['Bin'] == vocbin,'BUS_EXH'].values )
        if inputs.ld_diesel_on:
            ld_dsl_ef_df.loc[:,vocbin] = (ld_dsl_df.loc[:,'VOC_EXH']*voc_wt.loc[voc_wt['Bin'] == vocbin,'LD_DSL_EXH'].values/100/voc_mw.loc[voc_mw['Bin'] == vocbin,'LD_DSL_EXH'].values )
    
    gas_ef_df = gas_ef_df.reset_index()
    bus_ef_df = bus_ef_df.reset_index()
    dsl_ef_df = dsl_ef_df.reset_index()
    if inputs.ld_diesel_on:
        ld_dsl_ef_df = ld_dsl_ef_df.reset_index()
    gas_ef_df.loc[:,'Type'] = ['gas']*gas_ef_df.loc[:,'FIPS'].size
    dsl_ef_df.loc[:,'Type'] = ['bus']*dsl_ef_df.loc[:,'FIPS'].size
    bus_ef_df.loc[:,'Type'] = ['dsl']*bus_ef_df.loc[:,'FIPS'].size
    if inputs.ld_diesel_on:
        ld_dsl_ef_df.loc[:,'Type'] = ['ld_dsl']*ld_dsl_ef_df.loc[:,'FIPS'].size

    # Combine everything into a single dataframe
    if inputs.ld_diesel_on:
        emissionfactor_df = pd.concat([gas_ef_df.fillna(0),dsl_ef_df.fillna(0),bus_ef_df.fillna(0),ld_dsl_ef_df.fillna(0)])
    else:
        emissionfactor_df = pd.concat([gas_ef_df.fillna(0),dsl_ef_df.fillna(0),bus_ef_df.fillna(0)])
    
    return emissionfactor_df.reset_index(drop=True)
# end make_emissionfactors

# make_scalefactors
# makes scaling factors by state for each category
def make_scalefactors(year,month):
    ef_sf_fn_dict = inputs().get_input_filenames(year = year,month = month)
    ann_scale = pd.read_csv(inputs.input_dir+ef_sf_fn_dict['ann_fuelscale_fn'])
    month_fuel_scale = pd.read_csv(inputs.input_dir+ef_sf_fn_dict['monthly_fuel_fn'])
    if inputs.ld_diesel_on:
        month_fuel_scale.rename(mapper={'GAS_URB':'GAS_URB_FUEL','GAS_RUR':'GAS_RUR_FUEL',
                             'DSL_URB':'DSL_URB_FUEL','DSL_RUR':'DSL_RUR_FUEL',
                             'BUS_URB':'BUS_URB_FUEL','BUS_RUR':'BUS_RUR_FUEL',
                             'LD_DSL_URB':'LD_DSL_URB_FUEL','LD_DSL_RUR':'LD_DSL_RUR_FUEL',},inplace=True,axis='columns')
    else:
        month_fuel_scale.rename(mapper={'GAS_URB':'GAS_URB_FUEL','GAS_RUR':'GAS_RUR_FUEL',
                             'DSL_URB':'DSL_URB_FUEL','DSL_RUR':'DSL_RUR_FUEL',
                             'BUS_URB':'BUS_URB_FUEL','BUS_RUR':'BUS_RUR_FUEL'},inplace=True,axis='columns')
    
    month_scale = pd.read_csv(inputs.input_dir+inputs.monthly_fn)
    month_padd_scale = pd.read_csv(inputs.input_dir+ef_sf_fn_dict['monthly_PADDtoSTATE_fn'])
    
    month_scale_factors = month_scale[month_scale['Month']==month].reset_index(drop=True)
    df = month_padd_scale[month_padd_scale['Month']==month].reset_index(drop=True)
    df = df.join(month_fuel_scale[month_fuel_scale['Month']==month].reset_index(drop=True).set_index('PADD').drop(columns=['Month']),on='PADD')
    df = df.join(ann_scale.set_index('Name'),on='Name')

    df_urb = df[['Month','PADD','Name','FIPS']].copy().reset_index(drop=True)
    df_rur = df[['Month','PADD','Name','FIPS']].copy().reset_index(drop=True)
    
    #Construct full scaling factor
    df_urb['gas_sf'] = (df['state_scale']*df['GAS_URB_FUEL']*month_scale_factors['GAS_URB'].values*df['GAS']).fillna(1) 
    df_urb['dsl_sf'] = (df['state_scale']*df['DSL_URB_FUEL']*month_scale_factors['DSL_URB'].values*df['DSL']).fillna(1) 
    df_urb['bus_sf'] = (df['state_scale']*df['BUS_URB_FUEL']*month_scale_factors['BUS_URB'].values*df['BUS']).fillna(0)
    if inputs.ld_diesel_on:
        df_urb['ld_dsl_sf'] = (df['state_scale']*df['LD_DSL_URB_FUEL']*month_scale_factors['GAS_URB'].values*df['DSL']).fillna(1)
    df_urb['urb_rur'] = ['Urb']*df_urb.shape[0]
    
    df_rur['gas_sf'] = (df['state_scale']*df['GAS_RUR_FUEL']*month_scale_factors['GAS_RUR'].values*df['GAS']).fillna(1) 
    df_rur['dsl_sf'] = (df['DSL_RUR_FUEL']*month_scale_factors['DSL_RUR'].values*df['DSL']).fillna(1) 
    df_rur['bus_sf'] = (df['state_scale']*df['BUS_RUR_FUEL']*month_scale_factors['BUS_RUR'].values*df['BUS']).fillna(0)
    if inputs.ld_diesel_on:
        df_rur['ld_dsl_sf'] = (df['state_scale']*df['LD_DSL_RUR_FUEL']*month_scale_factors['GAS_RUR'].values*df['DSL']).fillna(1)
    df_rur['urb_rur'] = ['Rur']*df_rur.shape[0]
    
    df_complete = pd.concat([df_urb,df_rur])
    
    return df_complete.reset_index(drop=True)
# end make_scalefactors

# scale_fuel_ds
# takes fuel dataset and scales to the appropriate month
def scale_fuel_ds(fuel_ds,fuel_df,monthscale_df):
    
    fuel_ds_scaled = fuel_ds.copy(deep=True)
    # loop over the states available 
    for stfips in monthscale_df['FIPS'].drop_duplicates():
        if np.isnan(stfips):
            pass
        else:
            print('Scaling fuel for FIPS:'+ str(int(stfips)))
            # Find the indicies that correspond to urban and rural locations
            ind_x_urb = (fuel_df.loc[
                np.logical_and(fuel_df['state_fips']==stfips, fuel_df['urb_rur']=='Urb'),'col'
                               ].values.astype(int)).tolist()
            ind_y_urb = (fuel_df.loc[
                np.logical_and(fuel_df['state_fips']==stfips, fuel_df['urb_rur']=='Urb'),'row'
                               ].values.astype(int)).tolist()
            ind_x_rur = (fuel_df.loc[
                np.logical_and(fuel_df['state_fips']==stfips, fuel_df['urb_rur']=='Rur'),'col'
                               ].values.astype(int)).tolist()
            ind_y_rur = (fuel_df.loc[
                np.logical_and(fuel_df['state_fips']==stfips, fuel_df['urb_rur']=='Rur'),'row'
                               ].values.astype(int)).tolist()
            
            # Scale the fuel consumption for urban
            fuel_ds_scaled['gas'].data[ind_y_urb,ind_x_urb] = fuel_ds['gas'].data[ind_y_urb,ind_x_urb]*monthscale_df.loc[
                np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Urb'),'gas_sf'
                ].values
            fuel_ds_scaled['dsl'].data[ind_y_urb,ind_x_urb] = fuel_ds['dsl'].data[ind_y_urb,ind_x_urb]*monthscale_df.loc[
                np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Urb'),'dsl_sf'
                ].values
            fuel_ds_scaled['bus'].data[ind_y_urb,ind_x_urb] = fuel_ds['gas'].data[ind_y_urb,ind_x_urb]*monthscale_df.loc[
               np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Urb'),'bus_sf'
                ].values
            if inputs.ld_diesel_on:
                fuel_ds_scaled['ld_dsl'].data[ind_y_urb,ind_x_urb] = fuel_ds['ld_dsl'].data[ind_y_urb,ind_x_urb]*monthscale_df.loc[
                    np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Urb'),'ld_dsl_sf'
                    ].values
            
            # Scale the fuel consumption for rural
            fuel_ds_scaled['gas'].data[ind_y_rur,ind_x_rur] = fuel_ds['gas'].data[ind_y_rur,ind_x_rur]*monthscale_df.loc[
                np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Rur'),'gas_sf'
                ].values
            fuel_ds_scaled['dsl'].data[ind_y_rur,ind_x_rur] = fuel_ds['dsl'].data[ind_y_rur,ind_x_rur]*monthscale_df.loc[
                np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Rur'),'dsl_sf'
                ].values
            fuel_ds_scaled['bus'].data[ind_y_rur,ind_x_rur] = fuel_ds['gas'].data[ind_y_rur,ind_x_rur]*monthscale_df.loc[
               np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Rur'),'bus_sf'
                ].values
            if inputs.ld_diesel_on:
                fuel_ds_scaled['ld_dsl'].data[ind_y_rur,ind_x_rur] = fuel_ds['ld_dsl'].data[ind_y_rur,ind_x_rur]*monthscale_df.loc[
                    np.logical_and(monthscale_df['FIPS']==stfips,monthscale_df['urb_rur']=='Rur'),'ld_dsl_sf'
                    ].values
        # end if
    # # # end for loop
    
    return fuel_ds_scaled
# end scale_fuel_ds

# add_diurnal
def add_diurnal(fuel_ds,fuel_df,diurnal_df):
    
    # fuel_ds_full = fuel_ds.expand_dims({'day':np.arange(0,len(inputs.dow_list)),
    #                                'Time':np.arange(0,24)})
    
    extra_dims = xr.DataArray(data=np.zeros(shape=(len(inputs.dow_list),24)),
        dims=['day','Time'])
    
    fuel_ds_full = xr.broadcast(fuel_ds,extra_dims.to_dataset(name='extra_dims'))[0].copy(deep=True)
    # fuel_ds_full = fuel_ds_full.drop(['tz','urb_rur','fips'])
    fuel_ds_full = fuel_ds_full.drop(['fips'])
    
    # loop over the states available 
    for tz in diurnal_df['Timezone'].drop_duplicates():
        print('Adding diurnal profiles for Timezone:'+ str(tz))
        for day in inputs.dow_list:
            day_map = inputs.day_map
            # Find the indicies that correspond to urban and rural locations
            ind_x_urb = (fuel_df.loc[
                np.logical_and(fuel_df['tz']==tz, fuel_df['urb_rur']=='Urb'),'col'
                               ].values.astype(int)).tolist()
            ind_y_urb = (fuel_df.loc[
                np.logical_and(fuel_df['tz']==tz, fuel_df['urb_rur']=='Urb'),'row'
                               ].values.astype(int)).tolist()
            ind_x_rur = (fuel_df.loc[
                np.logical_and(fuel_df['tz']==tz, fuel_df['urb_rur']=='Rur'),'col'
                               ].values.astype(int)).tolist()
            ind_y_rur = (fuel_df.loc[
                np.logical_and(fuel_df['tz']==tz, fuel_df['urb_rur']=='Rur'),'row'
                               ].values.astype(int)).tolist()
            
            # Scale the fuel consumption for urban
            if ind_x_urb or ind_y_urb:
                da_gas_urb = xr.DataArray(data=diurnal_df.loc[
                    np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'GAS_URB']
                    .values,dims=['Time'])
                da_dsl_urb = xr.DataArray(data=diurnal_df.loc[
                    np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'DSL_URB']
                    .values,dims=['Time'])
                da_bus_urb = xr.DataArray(data=diurnal_df.loc[
                    np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'BUS_URB']
                    .values,dims=['Time'])
                if inputs.ld_diesel_on:
                    da_ld_dsl_urb = xr.DataArray(data=diurnal_df.loc[
                        np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'LD_DSL_URB']
                        .values,dims=['Time'])
                for hour in np.arange(0,da_gas_urb.size):
                    fuel_ds_full['gas'].data[ind_y_urb,ind_x_urb,day_map[day],hour] = fuel_ds['gas'].data[ind_y_urb,ind_x_urb]*da_gas_urb[hour].data/24
                    fuel_ds_full['dsl'].data[ind_y_urb,ind_x_urb,day_map[day],hour] = fuel_ds['dsl'].data[ind_y_urb,ind_x_urb]*da_dsl_urb[hour].data/24
                    fuel_ds_full['bus'].data[ind_y_urb,ind_x_urb,day_map[day],hour] = fuel_ds['bus'].data[ind_y_urb,ind_x_urb]*da_bus_urb[hour].data/24
                    if inputs.ld_diesel_on:
                        fuel_ds_full['ld_dsl'].data[ind_y_urb,ind_x_urb,day_map[day],hour] = fuel_ds['ld_dsl'].data[ind_y_urb,ind_x_urb]*da_ld_dsl_urb[hour].data/24
            else:
                pass
            
            # Scale the fuel consumption for rural
            if ind_x_rur or ind_y_rur:
                da_gas_rur = xr.DataArray(data=diurnal_df.loc[
                    np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'GAS_RUR']
                    .values,dims=['Time'])
                da_dsl_rur = xr.DataArray(data=diurnal_df.loc[
                    np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'DSL_RUR']
                    .values,dims=['Time'])
                da_bus_rur = xr.DataArray(data=diurnal_df.loc[
                    np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'BUS_RUR']
                    .values,dims=['Time'])
                if inputs.ld_diesel_on:
                    da_ld_dsl_rur = xr.DataArray(data=diurnal_df.loc[
                        np.logical_and(diurnal_df['Timezone']==tz, diurnal_df['Day']==day),'LD_DSL_RUR']
                        .values,dims=['Time'])
                for hour in np.arange(0,da_gas_rur.size):
                    fuel_ds_full['gas'].data[ind_y_rur,ind_x_rur,day_map[day],hour] = fuel_ds['gas'].data[ind_y_rur,ind_x_rur]*da_gas_rur[hour].data/24
                    fuel_ds_full['dsl'].data[ind_y_rur,ind_x_rur,day_map[day],hour] = fuel_ds['dsl'].data[ind_y_rur,ind_x_rur]*da_dsl_rur[hour].data/24
                    fuel_ds_full['bus'].data[ind_y_rur,ind_x_rur,day_map[day],hour] = fuel_ds['bus'].data[ind_y_rur,ind_x_rur]*da_bus_rur[hour].data/24
                    if inputs.ld_diesel_on:
                        fuel_ds_full['ld_dsl'].data[ind_y_rur,ind_x_rur,day_map[day],hour] = fuel_ds['ld_dsl'].data[ind_y_rur,ind_x_rur]*da_ld_dsl_rur[hour].data/24
            else:
                pass
            
        # end day loop
    # end for loop
    fuel_ds_full['gas'].attrs['units'] = 'kg/hr'
    fuel_ds_full['bus'].attrs['units'] = 'kg/hr'
    fuel_ds_full['dsl'].attrs['units'] = 'kg/hr'
    if inputs.ld_diesel_on:
        fuel_ds_full['ld_dsl'].attrs['units'] = 'kg/hr'
    

    if inputs.ld_diesel_on:
        fuels = ['gas','bus','dsl','ld_dsl']
    else:
        fuels = ['gas','bus','dsl']
    #for ii in fuels:
    #    fuel_ds_full[ii].encoding={'dtype': 'float32', 'chunksizes':(fuel_ds_full.sizes['south_north'], fuel_ds_full.sizes['west_east'], 1,1),
    #                      'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
    # # for ii in ['tz','fips','urb_rur']:
    # #     fuel_ds_full[ii].encoding={'chunksizes':(fuel_ds_full.sizes['south_north'], fuel_ds_full.sizes['west_east'], 1,1),
    # #                      'zlib': True, 'complevel': 1, '_FillValue': None}
    #fuel_ds_full[fuels].to_netcdf(inputs.output_base_dir+'fuel_full.nc',format='netCDF4',engine='netcdf4')
    return fuel_ds_full
# end add_diurnal

# make_full_emis
def write_full_emis(emissionfactor_df,fuel_ds_full,fuel_ds, year, month):
    
    year_dir_name = inputs.out_struct['year']['name'][year]
    # month_dir_name ='Month' + str(month).zfill(2)
    month_dir_name = inputs.out_struct['month']['name'][month]
    
    # fuel_types_list = emissionfactor_df['Type'].drop_duplicates() # different types of fuel
    fuel_types_list = inputs.output_fuel_types
    
    fips_list = emissionfactor_df['FIPS'].drop_duplicates()
    day_map = inputs.day_map
    fuel_ds_full = fuel_ds_full.transpose('day','Time','south_north','west_east') # transpose so will have same order as template
    
    # hours corresponding to each day half
    times_1= np.arange(0,12)    
    times_2= np.arange(12,24)
    
    # Species to output
    if inputs.all_species:
        species_list = xr.open_dataset(inputs.template_1).drop(labels='Times').data_vars
    else:
        species_list = inputs.species_sub_list
    
    # Loop over the days of interest
    for day in inputs.dow_list:
        # open different templates for different day halves
        if inputs.low_mem_on:
            template_1 = xr.open_dataset(inputs.template_1,chunks={'Time':1})
            template_2 = xr.open_dataset(inputs.template_2,chunks={'Time':1})
        else:
            template_1 = xr.load_dataset(inputs.template_1)
            template_2 = xr.load_dataset(inputs.template_2)
            # template_1 = xr.open_dataset(inputs.template_1, chunks = {'Time':1})
            # template_2 = xr.open_dataset(inputs.template_2, chunks = {'Time':1})
        # Loop over all species
        for species in species_list:
            print('Running species: ' + species + ' on day : ' + day)
            
            # Loop over the different fuel types
            for ind_fuel,fuel_type in enumerate(fuel_types_list):
                
                # if the emission factor is the same nationally, dont loop over all states
                if (inputs.dsl_ef_national and fuel_type in inputs.dsl_ef_national_type_list):
                    ef = emissionfactor_df.loc[np.logical_and(emissionfactor_df['Type']==fuel_type,emissionfactor_df['FIPS']==0),species].values
                    if ef == 0:
                        pass
                    da = fuel_ds_full[fuel_type][day_map[day],:,:,:]*ef
                else:
                    # if the emission factor varies by state, then loop over all states and then sum the result
                    for ind_fips,stfips in enumerate(fips_list):
                        state_mask = xr.where(fuel_ds['fips']==stfips,1,0)
                    
                        ef = emissionfactor_df.loc[np.logical_and(emissionfactor_df['Type']==fuel_type,emissionfactor_df['FIPS']==stfips),species].values
                        if ef == 0:
                            pass
                        if ind_fips == 0:
                            da = xr.where(state_mask,fuel_ds_full[fuel_type][day_map[day],:,:,:]*ef,0)
                        else:
                            da = da + xr.where(state_mask,fuel_ds_full[fuel_type][day_map[day],:,:,:]*ef,0)
                        # end if
                        state_mask.close()
                    #end for over states
                # end if for fuel types
                if ind_fuel ==0:
                    da_full = da
                else:
                    da_full = da_full + da
                del da
            # end for loop over fuel types
            da_full = da_full.transpose('Time','south_north','west_east')
            # if not an HC bin then convert grams to metric tons
            if species[:2] != 'HC':
                da_full.data = da_full.data/(10**6)
            # end for loop over states
            template_1[species].data = da_full[times_1,:,:].data
            template_2[species].data = da_full[times_2,:,:].data
            del da_full
            gc.collect()
        # end loop over species
        
        # Fill in time values
        date = str(year)+'-'+str(month).zfill(2)+'-'+'15_'
        hrs_00z = ['00:00:00','01:00:00','02:00:00','03:00:00','04:00:00','05:00:00',
                   '06:00:00','07:00:00','08:00:00','09:00:00','10:00:00','11:00:00']
        hrs_12z = ['12:00:00','13:00:00','14:00:00','15:00:00','16:00:00','17:00:00',
                   '18:00:00','19:00:00','20:00:00','21:00:00','22:00:00','23:00:00']
        dates_00z = [i+j for i,j in zip([date]*12,hrs_00z)]
        dates_12z = [i+j for i,j in zip([date]*12,hrs_12z)]
        template_1['Times'].data = dates_00z
        template_2['Times'].data = dates_12z
        
        #If file output turned on, write to disk with the following encoding options
        if inputs.output_on:
            for ii in species_list:
                template_1[ii].encoding={'dtype': 'float32', 'chunksizes':(1,fuel_ds_full.sizes['south_north'], fuel_ds_full.sizes['west_east']),
                          'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
                template_2[ii].encoding={'dtype': 'float32', 'chunksizes':(1,fuel_ds_full.sizes['south_north'], fuel_ds_full.sizes['west_east']),
                          'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
                # template_1[ii].encoding=inputs.out_file_species_encoding
                # template_2[ii].encoding=inputs.out_file_species_encoding
            
            template_1['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}
            template_1['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}
            template_1['Times'].encoding={'char_dim_name':'DateStrLen'}
            template_2['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}
            template_2['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}
            template_2['Times'].encoding={'char_dim_name':'DateStrLen'}
            
            out_dir = inputs.output_base_dir  + year_dir_name + '/' + month_dir_name + '/' + day + '/'
    
            out_fn_1 = out_dir + inputs.out_fn1#'onroad_00to12Z.nc'
            out_fn_2 = out_dir + inputs.out_fn2#'onroad_12to24Z.nc'
            print('Writing file: ', out_fn_1)
            template_1.to_netcdf(out_fn_1,format='netCDF4',engine='netcdf4')
            del template_1
            print('Writing file: ', out_fn_2)
            template_2.to_netcdf(out_fn_2,format='netCDF4',engine='netcdf4')
            del template_2
            gc.collect()
            
    # end loop over days
        
# end make_full_emis

def main():
    pd.options.mode.chained_assignment = 'raise'
    start_time = time.time()
    print('START OF PROGRAM: ')
    print('--------------------------------------------------')
    if inputs.output_on:
        print('Output files: True')
        print('Selected fuel types for output:' + str(inputs.output_fuel_types))
    else:
        print('Output files: False')
    if inputs.all_species:
        print('Output Species: All available species')
    else:
        print('Output Species: ' + str(inputs.species_sub_list))
    
    
    # get unscaled fuel dataset from shape file
    fuel_ds, fuel_df = initialize_fuel_ds() 
    for y in inputs.years:
        for m in inputs.months:
            print('Calculating emissions for Year: ' + str(y) + ',  Month: ' + str(m))
            
            start_month_time = time.time()
            # Get diurnal scalings, monthly scalings, and emission factors for the time period
            diurnalscale_df = make_diurnalfactors(month=m)
            emissionfactor_df = make_emissionfactors(year = y,month = m)
            
            monthscale_df = make_scalefactors(year = y,month = m)
            # save overall scalings to file if option on
            if inputs.save_fuel_scaling:
                monthscale_df.to_csv(inputs.output_base_dir  + inputs.out_struct['year']['name'][y] + '/' + inputs.out_struct['month']['name'][m] + '/' + inputs.fuel_scaling_fn,index=False)
            
            # Scale the fuel consumption for each category to the appropriate time period
            fuel_ds_scaled  = scale_fuel_ds(fuel_ds = fuel_ds,fuel_df = fuel_df,monthscale_df = monthscale_df)
            # save to file if option on
            if inputs.save_monthscaled_fuel:
                for ii in ['gas','bus','dsl']:
                    fuel_ds_scaled[ii].encoding={'dtype': 'float32', 'chunksizes':(fuel_ds.sizes['south_north'], fuel_ds.sizes['west_east']),
                                      'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
                fuel_ds_scaled['fips'].encoding={'chunksizes':(fuel_ds.sizes['south_north'], fuel_ds.sizes['west_east']),
                                  'zlib': True, 'complevel': 1, '_FillValue': None}
                fuel_ds_scaled.to_netcdf(inputs.output_base_dir  + inputs.out_struct['year']['name'][y] + '/' + 
                                         inputs.out_struct['month']['name'][m] + '/' + inputs.monthscaled_fuel_fn,format='netCDF4',engine='netcdf4')
            
            # Add diurnal variability for each timezone
            fuel_ds_full = add_diurnal(fuel_ds=fuel_ds_scaled,fuel_df=fuel_df,diurnal_df=diurnalscale_df )
            # Save to file if option on
            if inputs.save_diurnal_fuel:
                for ii in ['gas','bus','dsl']:
                    fuel_ds_full[ii].encoding={'dtype': 'float32', 'chunksizes':(1,1,fuel_ds.sizes['south_north'], fuel_ds.sizes['west_east']),
                                      'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
                #fuel_ds_full['fips'].encoding={'chunksizes':(fuel_ds.sizes['south_north'], fuel_ds.sizes['west_east']),
                #                  'zlib': True, 'complevel': 1, '_FillValue': None}
                fuel_ds_full.to_netcdf(inputs.output_base_dir  + inputs.out_struct['year']['name'][y] + '/' + 
                                         inputs.out_struct['month']['name'][m] + '/' + inputs.diurnal_fuel_fn,format='netCDF4',engine='netcdf4')
            
            
            # Applies emission factors and writes the relevant files to disk
            write_full_emis(emissionfactor_df=emissionfactor_df,fuel_ds_full=fuel_ds_full,fuel_ds=fuel_ds,
                           year=y,month=m)
            
            # Delete the variables that wont be used in the next iteration 
            del fuel_ds_scaled, fuel_ds_full, diurnalscale_df, emissionfactor_df, monthscale_df
            gc.collect()
            end_month_time = time.time()
            print('Time to make files for Year : ' + str(y) + ',  Month: ' + str(m) + ' is ' + str( round(end_month_time-start_month_time ) )  + ' s')
            
    # print(fuel_ds)
    print('END OF PROGRAM: ')
    print('--------------------------------------------------')
    end_time = time.time()
    print('Total elapsed time (s): ' + str(round(end_time-start_time)))
    
# End Main

if __name__ == "__main__":
    main()

    



    

