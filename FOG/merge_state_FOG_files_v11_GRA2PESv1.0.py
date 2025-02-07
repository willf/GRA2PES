# Uses geopandas_env
import numpy as np
import pandas as pd
import geopandas as gpd
import time
import xarray as xr


class inputs():
    #input_dir = '/wrk/d2/charkins/FOG_data/states_output/'
    #output_dir = '/wrk/d2/charkins/FOG_data/states_output/'
    #ncf_output_dir = '/wrk/d2/charkins/FOG_data/FOG_only_ncf/'
    
    input_dir = '/wrk/users/charkins/emissions/FOG/FOG_data/states_output_v11_GRA2PESv1.0/'
    output_dir = '/wrk/users/charkins/emissions/FOG/FOG_data/states_output_v11_GRA2PESv1.0/'
    ncf_output_dir = '/wrk/users/charkins/emissions/FOG/FOG_data/FOG_ncf_v11_GRA2PESv1.0/'
    
    # state_list = ['CO','WA']
    # state_list = ['WA','CO']

    
    state_list = ['AL','AR','AZ','CA','CO','FL',
                  'KS','KY','LA','MD','MI','MN',
                  'MO','MS','MT','ND','NE','NM',
                  'NV','NY','OH','OK','OR','PA',
                  'SD','TN','TX','UT','VA','WA',
                  'WV','WY']
    
    combine_cols = ['Oprod','Gprod','DrCount','SUM_HtdryC','SUM_HtwetC','Sum_WcCO2','SUM_LfCO2','LeasEngCO2','DyFCO2','LfFCO2',
                    'HtFCO2','LcFCO2','WcFCO2','DrFCO2','TotFCO2','LcFNOx','WcFNOx','LfFNOx','HtFNOx','DyFNOx',
                    'DrFNOx','TotFNOx','flare_NOx','NoppFNOx','CH4','ALK1','ALK4','ALK5','ARO1','BENZENE','BUTANES',
                    'PENTANES','TOLUENE','M_XYLENE','O_XYLENE','P_XYLENE','PROPANE','NMHC',
                    'CO','NH3','PM10-PRI','PM25-PRI','SO2','ETHENE','HCHO','CCHO','ACETYLENE',
                    'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                    'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']
    
    remap_names = {'TotFCO2':'CO2',
                   'NoppFNOx':'NOX',
                   'CH4':'HC01',
                   'ALK1':'HC02',
                   'ALK4':'HC05',
                   'ALK5':'HC06',
                   'ETHENE':'HC07',
                   'ARO1':'HC12',
                   'HCHO':'HC14',
                   'CCHO':'HC15',
                   'ACETYLENE':'HC37',
                   'BENZENE':'HC38',
                   'BUTANES':'HC39',
                   'PENTANES':'HC40',
                   'TOLUENE':'HC41',
                   'M_XYLENE':'HC42',
                   'O_XYLENE':'HC43',
                   'P_XYLENE':'HC44',
                   'PROPANE':'HC45',
                   }
    # variables to add to netcdf template
    ncf_vars = list(remap_names.values()) + ['CO','NH3','PM10-PRI','PM25-PRI','SO2',
                    'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                    'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19','VOC']
    
    in_fn_base = '[ST]_final_data_YYYY_Month##.'
    out_fn_base = 'FOG_data_YYYY_Month##.'
    ncf_out_fn_base = 'FOG'
    
    # Domain info
    domain_fn='/wrk/charkins/Domains/'+'nei04k_domain.shp'
    
    # Netcdf Template locations
    template_fn = ['/wrk/charkins/FIVE/FIVE_onroad_python/template/onroad_00to12Z.nc', '/wrk/charkins/FIVE/FIVE_onroad_python/template/onroad_12to24Z.nc' ]
    
    # Molecular weight info
    voc_convert_list = ['HC01','HC02','HC05','HC06','HC07','HC12','HC14','HC15','HC37','HC38','HC39','HC40','HC41','HC42','HC43','HC44','HC45']
    mw_fn = '/wrk/charkins/emissions/FOG/inputs/FOG_mw.csv'
    
    years=[2023]
    months=[1,2,3,4,5,6,7,8,9,10,11,12]
    
    #options
    save_gpkg = False
    save_pkl = True
    make_ncf = True
    save_ncf = True
    save_state_totals = True
    state_totals_dir = '/wrk/charkins/emissions/FOG/merge_states_summary/'
    state_totals_fn = 'FOG_state_totals_2021_v11_GRA2PESv1.0.csv'
    
# end inputs class

# adds state data to the full domain
def add_state_data(full_df,state_df,state):
    full_df=full_df.copy(deep=True)
    state_df=state_df.copy(deep=True)
    
    # Set the index to ij_key so that we can get the info we want
    full_df['Ind'] = (full_df['Row']*(full_df['Col'].max()+1))+(full_df['Col'])
    state_df['Ind'] = (state_df['Row']*(full_df['Col'].max()+1))+(state_df['Col'])
    full_df = full_df.set_index('Ind')
    state_df = state_df.set_index('Ind')[inputs.combine_cols]  # select only the columns we want to combine

    for col in inputs.combine_cols:
        full_df.loc[state_df.index,col] = full_df.loc[state_df.index,col].add(state_df.loc[:,col],fill_value=0) 

    full_df.reset_index(inplace=True,drop=True)
    
    return full_df
# end add_state_data


def add_columns(df):
    df = df.copy(deep=True)
    
    df['Oprod']=np.zeros(df['ij_key'].shape)
    df['Gprod']=np.zeros(df['ij_key'].shape)
    df['DrCount']=np.zeros(df['ij_key'].shape)
    df['SUM_HtdryC']=np.zeros(df['ij_key'].shape)
    df['SUM_HtwetC']=np.zeros(df['ij_key'].shape)
    df['Sum_WcCO2']=np.zeros(df['ij_key'].shape)
    df['SUM_LfCO2']=np.zeros(df['ij_key'].shape)
    df['LeasEngCO2']=np.zeros(df['ij_key'].shape)
    df['DyFCO2']=np.zeros(df['ij_key'].shape)
    df['LfFCO2']=np.zeros(df['ij_key'].shape)
    df['HtFCO2']=np.zeros(df['ij_key'].shape)
    df['LcFCO2']=np.zeros(df['ij_key'].shape)
    df['WcFCO2']=np.zeros(df['ij_key'].shape)
    df['DrFCO2']=np.zeros(df['ij_key'].shape)
    df['TotFCO2']=np.zeros(df['ij_key'].shape)
    df['LcFNOx']=np.zeros(df['ij_key'].shape)
    df['flare_NOx']=np.zeros(df['ij_key'].shape)
    df['WcFNOx']=np.zeros(df['ij_key'].shape)
    df['LfFNOx']=np.zeros(df['ij_key'].shape)
    df['HtFNOx']=np.zeros(df['ij_key'].shape)
    df['DyFNOx']=np.zeros(df['ij_key'].shape)
    df['DrFNOx']=np.zeros(df['ij_key'].shape)
    df['TotFNOx']=np.zeros(df['ij_key'].shape)
    df['NoppFNOx']=np.zeros(df['ij_key'].shape)
    
    # Add all VOC species as zero values
    df['CH4']=np.zeros(df['ij_key'].shape)
    df['ALK1']=np.zeros(df['ij_key'].shape)
    df['ALK4']=np.zeros(df['ij_key'].shape)
    df['ALK5']=np.zeros(df['ij_key'].shape)
    df['ARO1']=np.zeros(df['ij_key'].shape)
    df['BENZENE']=np.zeros(df['ij_key'].shape)
    df['BUTANES']=np.zeros(df['ij_key'].shape)
    df['PENTANES']=np.zeros(df['ij_key'].shape)
    df['TOLUENE']=np.zeros(df['ij_key'].shape)
    df['M_XYLENE']=np.zeros(df['ij_key'].shape)
    df['O_XYLENE']=np.zeros(df['ij_key'].shape)
    df['P_XYLENE']=np.zeros(df['ij_key'].shape)
    df['PROPANE']=np.zeros(df['ij_key'].shape)
    df['NMHC']=np.zeros(df['ij_key'].shape)

    # Add all the rest of the species added based on ratio to CO2
    df['CO']=np.zeros(df['ij_key'].shape)
    df['NH3']=np.zeros(df['ij_key'].shape)
    df['PM10-PRI']=np.zeros(df['ij_key'].shape)
    df['PM25-PRI']=np.zeros(df['ij_key'].shape)
    df['SO2']=np.zeros(df['ij_key'].shape)
    df['ETHENE']=np.zeros(df['ij_key'].shape) # also HC07
    df['HCHO']=np.zeros(df['ij_key'].shape)# also HC14
    df['CCHO']=np.zeros(df['ij_key'].shape)# also HC15
    df['ACETYLENE']=np.zeros(df['ij_key'].shape)# also HC37
    df['PM01']=np.zeros(df['ij_key'].shape)
    df['PM02']=np.zeros(df['ij_key'].shape)
    df['PM03']=np.zeros(df['ij_key'].shape)
    df['PM04']=np.zeros(df['ij_key'].shape)
    df['PM05']=np.zeros(df['ij_key'].shape)
    df['PM06']=np.zeros(df['ij_key'].shape)
    df['PM07']=np.zeros(df['ij_key'].shape)
    df['PM08']=np.zeros(df['ij_key'].shape)
    df['PM09']=np.zeros(df['ij_key'].shape)
    df['PM10']=np.zeros(df['ij_key'].shape)
    df['PM11']=np.zeros(df['ij_key'].shape)
    df['PM12']=np.zeros(df['ij_key'].shape)
    df['PM13']=np.zeros(df['ij_key'].shape)
    df['PM14']=np.zeros(df['ij_key'].shape)
    df['PM15']=np.zeros(df['ij_key'].shape)
    df['PM16']=np.zeros(df['ij_key'].shape)
    df['PM17']=np.zeros(df['ij_key'].shape)
    df['PM18']=np.zeros(df['ij_key'].shape)
    df['PM19']=np.zeros(df['ij_key'].shape)

    return df
    
# end add_columns

def save_output(inputs, gdf, year,month):
    
    if inputs.save_gpkg:
        out_fn = inputs.output_dir + inputs.out_fn_base.replace('YYYY',str(year)).replace('##',str(month).zfill(2)) + 'gpkg'
        print('Writing output gpkg for [Year: '+ str(year) + ', Month: '+ str(month)+ ']')
        gdf.to_file(out_fn,driver="GPKG")
    if inputs.save_pkl:
        out_fn = inputs.output_dir + inputs.out_fn_base.replace('YYYY',str(year)).replace('##',str(month).zfill(2))+ 'pkl.gz'
        print('Writing output pkl for [Year: '+ str(year) + ', Month: '+ str(month)+ ']')
        compression_opts={'method': 'gzip', 'compresslevel': 1}
        gdf.drop(columns='geometry').to_pickle(out_fn,compression=compression_opts)
        
# # end save_output

# Make FOG netcdf file to be saved out
def make_ncf(inputs,df,year,month):
    # ind_x = (df['Col'].values.astype(int)).tolist()
    # ind_y = (df['Row'].values.astype(int)).tolist()
    # print('Calc inds')
    ind_x = xr.DataArray((df['Col'].values.astype(int)).tolist(),dims=["d"])
    ind_y = xr.DataArray((df['Row'].values.astype(int)).tolist(),dims=["d"])
    
    # Loop over the two day half templates
    for template in inputs.template_fn:
        ds = xr.open_dataset(template)
        
        # Loop over the variables for which we have FOG data
        for var in inputs.ncf_vars:
            print('Assigning for var: ' + var )
            for t in np.arange(0,12):
                # convert variables from mol/day to mol/hr and mt/d to mt/hr
                ds[var].data[t,ind_y,ind_x] = df[var].astype(float).values/24
        
        # Add correct date strings for times variable
        date = str(year)+'-'+str(month).zfill(2)+'-'+'15_'
        
        if '00to12Z' in template: 
            hrs = ['00:00:00','01:00:00','02:00:00','03:00:00','04:00:00','05:00:00',
                       '06:00:00','07:00:00','08:00:00','09:00:00','10:00:00','11:00:00']
            half = 1
        elif '12to24Z' in template:
            hrs = ['12:00:00','13:00:00','14:00:00','15:00:00','16:00:00','17:00:00',
                       '18:00:00','19:00:00','20:00:00','21:00:00','22:00:00','23:00:00']
            half=2
            
        dates = [i+j for i,j in zip([date]*12,hrs)]
        ds['Times'].data = dates
        
        if inputs.save_ncf:
            save_ncf(inputs,ds,half,year,month)
            
# end make_ncf

# save out the FOG netcdf file 
def save_ncf(inputs,ds,half,year,month):
    
    import shutil
    
    # Set all encoding settings properly
    encoding_dict = {'dtype': 'float32', 'chunksizes':(1,ds.sizes['south_north'], ds.sizes['west_east']),
              'zlib': True, 'complevel': 1, '_FillValue': None , 'coordinates':'XLONG XLAT'}
    
    ncf_vars_all = ds.drop(labels='Times').data_vars
    for var in ncf_vars_all:
        ds[var].encoding=encoding_dict
    
    ds['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}
    ds['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}
    ds['Times'].encoding={'char_dim_name':'DateStrLen'}
    
    if half == 1: 
        #out_fn = inputs.ncf_output_dir + str(year) + '/Month' + str(month).zfill(2) + '/' + 'weekdy/' + inputs.ncf_out_fn_base.replace('YYYY',str(year)).replace('##',str(month).zfill(2)) + '_00to12Z.nc'
        out_fn = inputs.ncf_output_dir + str(year) + '/Month' + str(month).zfill(2) + '/' + 'weekdy/' + inputs.ncf_out_fn_base + '_00to12Z.nc'
    elif half == 2:
        #out_fn = inputs.ncf_output_dir + str(year) + '/Month' + str(month).zfill(2) + '/' + 'weekdy/' + inputs.ncf_out_fn_base.replace('YYYY',str(year)).replace('##',str(month).zfill(2)) + '_12to24Z.nc'
        out_fn = inputs.ncf_output_dir + str(year) + '/Month' + str(month).zfill(2) + '/' + 'weekdy/' + inputs.ncf_out_fn_base + '_12to24Z.nc'
    
    print('Writing file: ', out_fn)
    ds.to_netcdf(out_fn,format='netCDF4',engine='netcdf4')
    
    # Copy files for saturday and sunday
    shutil.copy(out_fn,out_fn.replace('weekdy/','satdy/'))
    shutil.copy(out_fn,out_fn.replace('weekdy/','sundy/'))
    
# end save_ncf

# converts variables from mol/time to mt/time
def convert_units_mt(inputs,df,var_list):
    df_out = df.copy(deep=True)
    mw_list = pd.read_csv(inputs.mw_fn)
    #Convert VOC units from mol/time to mt/time
    for var in var_list:
        mw = mw_list.loc[mw_list['Bin']==var,'OnG'].values[0]
        df_out[var]=df_out[var]*mw/1000000
    return df_out
# end convert units mt



def main():
    start_time = time.time()
    print('START OF PROGRAM: ')
    print('--------------------------------------------------')

    warning = None
    print('Reading in national domain....')
    
    # Read in shapefiles for domain so it only happens once
    domain_full=gpd.read_file(inputs.domain_fn)
    domain_full['ij_key']=domain_full['Row'].astype(str)+'.'+domain_full['Col'].astype(str)
    
    # Add all relevant columns to file
    domain_full = add_columns(domain_full)
    
    if inputs.save_state_totals:
        state_totals = pd.DataFrame({'Year' : [],
                               'Month' : [],
                               'State' : [],
                               'CO2':[],
                               'NOX':[],
                               'HC01':[],
                               'HC02':[],
                               'HC05':[],
                               'HC06':[],
                               'HC07':[],
                               'HC12':[],
                               'HC14':[],
                               'HC15':[],
                               'HC37':[],
                               'HC38':[],
                               'HC39':[],
                               'HC40':[],
                               'HC41':[],
                               'HC42':[],
                               'HC43':[],
                               'HC44':[],
                               'HC45':[],
                               'WcFNOx':[],
                               'LcFNOx':[],
                               'flare_NOx':[],
                               'DyFNOx':[],
                               'DrFNOx':[],
                               'HtFNOx':[],
                               'LfFNOx':[],
                               'DyFCO2':[],
                               'LfFCO2':[],
                               'HtFCO2':[],
                               'LcFCO2':[],
                               'WcFCO2':[],
                               'DrFCO2':[],
                               'CO':[],
                               'NH3':[],
                               'PM10-PRI':[],
                               'PM25-PRI':[],
                               'SO2':[],
                               'PM01':[],
                               'PM02':[],
                               'PM03':[],
                               'PM04':[],
                               'PM05':[],
                               'PM06':[],
                               'PM07':[],
                               'PM08':[],
                               'PM09':[],
                               'PM10':[],
                               'PM11':[],
                               'PM12':[],
                               'PM13':[],
                               'PM14':[],
                               'PM15':[],
                               'PM16':[],
                               'PM17':[],
                               'PM18':[],
                               'PM19':[],
                               })
    # end if
    
    # Loop over years in list
    for year in inputs.years:
        # Loop over months in list
        for month in inputs.months: 
            start_month_time = time.time()
            print('Working on Emissions for: ' + str(year) + ', Month : ' + str(month) )
            
            FOG_gdf = domain_full.copy(deep=True)
            
            # Loop over states in list
            for state in inputs.state_list:
                print('Adding State emissions files for: ' +state)
                state_df_fn = inputs.output_dir + inputs.in_fn_base.replace('[ST]',state).replace('YYYY',str(year)).replace('##',str(month).zfill(2))+ 'pkl.gz'
                state_df = pd.read_pickle(state_df_fn,compression={'method':'gzip'})
                
                if inputs.save_state_totals:
                    totals_append = pd.DataFrame({'Year' : [year],
                                           'Month' : [month],
                                           'State' : [state],
                                           'CO2':[state_df['TotFCO2'].sum()], #mt/d
                                           'NOX':[state_df['NoppFNOx'].sum()],#mt/d
                                           'HC01':[state_df['CH4'].sum()],  #mol/d
                                           'HC02':[state_df['ALK1'].sum()], #mol/d
                                           'HC05':[state_df['ALK4'].sum()], #mol/d
                                           'HC06':[state_df['ALK5'].sum()], #mol/d
                                           'HC07':[state_df['ETHENE'].sum()], #mol/d
                                           'HC12':[state_df['ARO1'].sum()], #mol/d
                                           'HC14':[state_df['HCHO'].sum()], #mol/d
                                           'HC15':[state_df['CCHO'].sum()], #mol/d
                                           'HC37':[state_df['ACETYLENE'].sum()], #mol/d
                                           'HC38':[state_df['BENZENE'].sum()], #mol/d
                                           'HC39':[state_df['BUTANES'].sum()], #mol/d
                                           'HC40':[state_df['PENTANES'].sum()], #mol/d
                                           'HC41':[state_df['TOLUENE'].sum()], #mol/d
                                           'HC42':[state_df['M_XYLENE'].sum()], #mol/d
                                           'HC43':[state_df['O_XYLENE'].sum()], #mol/d
                                           'HC44':[state_df['P_XYLENE'].sum()], #mol/d
                                           'HC45':[state_df['PROPANE'].sum()], #mol/d
                                           'WcFNOx':[state_df['WcFNOx'].sum()], #mt/d
                                           'LcFNOx':[state_df['LcFNOx'].sum()], #mt/d
                                           'flare_NOx':[state_df['flare_NOx'].sum()], # mt/d
                                           'DyFNOx':[state_df['DyFNOx'].sum()],#mt/d
                                           'DrFNOx':[state_df['DrFNOx'].sum()],#mt/d
                                           'HtFNOx':[state_df['HtFNOx'].sum()], #mt/d
                                           'LfFNOx':[state_df['LfFNOx'].sum()],#mt/d
                                           'DyFCO2':[state_df['DyFCO2'].sum()], #mt/d
                                           'LfFCO2':[state_df['LfFCO2'].sum()], #mt/d
                                           'HtFCO2':[state_df['HtFCO2'].sum()], #mt/d
                                           'LcFCO2':[state_df['LcFCO2'].sum()], #mt/d
                                           'WcFCO2':[state_df['WcFCO2'].sum()], #mt/d
                                           'DrFCO2':[state_df['DrFCO2'].sum()], #mt/d
                                           'CO':[state_df['CO'].sum()], #mt/d
                                           'NH3':[state_df['NH3'].sum()], #mt/d
                                           'PM10-PRI':[state_df['PM10-PRI'].sum()], #mt/d
                                           'PM25-PRI':[state_df['PM25-PRI'].sum()], #mt/d
                                           'SO2':[state_df['SO2'].sum()], #mt/d
                                           'PM01':[state_df['PM01'].sum()], #mt/d
                                           'PM02':[state_df['PM02'].sum()], #mt/d
                                           'PM03':[state_df['PM03'].sum()], #mt/d
                                           'PM04':[state_df['PM04'].sum()], #mt/d
                                           'PM05':[state_df['PM05'].sum()], #mt/d
                                           'PM06':[state_df['PM06'].sum()], #mt/d
                                           'PM07':[state_df['PM07'].sum()], #mt/d
                                           'PM08':[state_df['PM08'].sum()], #mt/d
                                           'PM09':[state_df['PM09'].sum()], #mt/d
                                           'PM10':[state_df['PM10'].sum()], #mt/d
                                           'PM11':[state_df['PM11'].sum()], #mt/d
                                           'PM12':[state_df['PM12'].sum()], #mt/d
                                           'PM13':[state_df['PM13'].sum()], #mt/d
                                           'PM14':[state_df['PM14'].sum()], #mt/d
                                           'PM15':[state_df['PM15'].sum()], #mt/d
                                           'PM16':[state_df['PM16'].sum()], #mt/d
                                           'PM17':[state_df['PM17'].sum()], #mt/d
                                           'PM18':[state_df['PM18'].sum()], #mt/d
                                           'PM19':[state_df['PM19'].sum()], #mt/d
                                           })
                    state_totals = pd.concat([state_totals,totals_append])
                
                FOG_gdf = add_state_data(FOG_gdf,state_df,state)
                


                # warning = write_warning_text(warning_eia, warning_flare, warning_voc, warning)

            # end for loop over states
            
            FOG_gdf.rename(columns=inputs.remap_names,inplace=True)
            
            # Convert Units to metric tons for VOC total
            FOG_gdf_mt = convert_units_mt(inputs,FOG_gdf,inputs.voc_convert_list)
            
            # Add total mass VOC
            FOG_gdf['VOC'] = 0
            for var in inputs.voc_convert_list:
                if var != 'HC01':
                    FOG_gdf['VOC'] = FOG_gdf['VOC']+FOG_gdf_mt[var]
                
            # Save the full output as gpkg and/or pickle
            save_output(inputs, FOG_gdf,year,month)
            # make netcdf and save the variable data to netcdf
            if inputs.make_ncf:
                make_ncf(inputs,FOG_gdf,year,month)
            
                
            end_month_time = time.time()
            print('Time elapsed to make file  for Year : ' + str(year) + ', Month : ' + str(month) + ', State: ' + state + ' is ' + str( round(end_month_time-start_month_time ) )  + ' s')
            # end Loop over months
        # end Loop over years
    # end for loop over state
    state_totals_mt = convert_units_mt(inputs,state_totals,inputs.voc_convert_list)
    state_totals['VOC'] = 0
    for var in inputs.voc_convert_list:
        if var != 'HC01':
            state_totals['VOC'] = state_totals['VOC']+state_totals_mt[var]
        
    state_totals.to_csv(inputs.state_totals_dir + inputs.state_totals_fn,index=False)
    
    print('END OF PROGRAM: ')
    print('--------------------------------------------------')
    end_time = time.time()
    print('Total elapsed time (s): ' + str(round(end_time-start_time)))
    
    # If there are warnings print them
    if (warning is not None):
        print('WARNINGS: ')
        print('--------------------------------------------------')
        print(warning['Message'])
        print('--------------------------------------------------')
    
    #print('NOTES: ')

    

if __name__ == "__main__":
    main()
