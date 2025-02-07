# V9 add columns to track drill specific CH4 and NMVOCs
# V10 adds complete CO2 
# V11 adds ratios to CO2 for CO, NH3, PM10-PRI, PM25-PRI, SO2, HC07 (ethylene), HC14 (formaldehyde), HC15 (acetaldehyde), HC37 (Acetylene)
# V11 also adds CH4 intensity trends over time

# Uses geopandas_env
import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
import time
import calendar
from os.path import exists
from os.path import basename
import glob

class inputs():
    input_dir = '/wrk/charkins/emissions/FOG/inputs/'
    data_dir = '/wrk/users/charkins/emissions/FOG/FOG_data/FOG_inputs/'
    output_dir = '/wrk/users/charkins/emissions/FOG/FOG_data/states_output_v11_GRA2PESv1.0/'
    
    #state_list = ['CO']
    
    state_list = ['AL','AR','AZ','CA','CO','FL',
                    'KS','KY','LA','MD','MI','MN',
                    'MO','MS','MT','ND','NE','NM',
                    'NV','NY','OH','OK','OR','PA',
                    'SD','TN','TX','UT','VA','WA',
                    'WV','WY']
    
    shape_dir = '/wrk/users/charkins/emissions/FOG/FOG_data/FOG_inputs/Every_map/'
    
    well_data_fn_base =  data_dir+'Production/Production_formatted/'+'[ST]_Monthly_Production_YYYY_Month##.csv'
    prod_summary_fn_base =  data_dir+'Production/Production_formatted/'+'[ST]_Production_totals_YYYY.csv'
    drill_data_fn_base =  data_dir+'Drilling/'+'[ST] Wells Table.CSV'
    flare_data_fn_base =  data_dir+'Flaring/'+'raw_flare_YYYY.csv'
    spatial_info_fn = input_dir + 'Spatial_info.xlsx'
    eia_info_fn = input_dir + 'eia.xlsx'
    basin_info_fn_base = input_dir + 'Basin_Ratio_CH20241212.xlsx'
    state_ratios_fn = data_dir + '2015_StateRatios.csv'
    
    out_fn_base = '[ST]_final_data_YYYY_Month##.'
    
    # Domain info
    # state_domains_dir: Directory to save state specific domain files if they dont exist
    state_domains_dir = shape_dir + 'state_domains/' 
    
    domain_fn='/wrk/charkins/Domains/'+'nei04k_domain.shp'
    basins_fn=shape_dir+'US_basins.shp'
    counties_fn=shape_dir+'USA_Counties.shp'
    states_fn=shape_dir+'USA_States_Generalized.shp'
    plays_fn='/wrk/charkins/emissions/FOG/inputs/methane_trends/US_plays_for_methane_trends.shp'
    play_methane_trends_fn='/wrk/charkins/emissions/FOG/inputs/methane_trends/US_plays_methane_trends_McDonald_EM2023_polynomial_fit.xlsx'
    
    state_fips_table_fn = input_dir+'us-state-ansi-fips.csv'
    # state_fips_table = pd.read_csv(state_fips_table_fn,
    #                                dtype={'stname':str,'st':str,'stusps':str},
    #                                skipinitialspace=True)
    state_fips_table = pd.read_csv(state_fips_table_fn,
                                   sep=',',
                                   dtype={'stname':str,'st':str,'stusps':str},
                                   skipinitialspace=True)
    
    years=[2023]
    months=[1,2,3,4,5,6,7,8,9,10,11,12]
    
    #options
    overwrite_state_domains = False # If true save new state domains regardless of whether they exist
    save_gpkg = False
    save_pkl = True
    
# end inputs class



# Reads domain files and cuts them down to be specific to the state
def make_state_domain(inputs,state,domain_full,counties_full,states_full):
    from os.path import exists
    
    print('Making State-specific domain files for: '+ state)
    
    st_fips = inputs.state_fips_table.loc[inputs.state_fips_table['stusps']==state,'st'].values[0]
    state_sub = states_full.loc[states_full['STATE_ABBR']==state,:]
    counties_sub=counties_full.loc[counties_full['STATE_FIPS']==st_fips,:]
    
    # If domain file specific to state doesn't already exist then create it 
    state_domain_fn = inputs.state_domains_dir +state+'_domain.gpkg'
    if (exists(state_domain_fn) and not inputs.overwrite_state_domains):
        domain_sub=gpd.read_file(state_domain_fn)
    else:
        domain_full=domain_full.copy(deep=True)
        domain_full['ij_key']=domain_full['Row'].astype(str)+'.'+domain_full['Col'].astype(str)
        domain_centpoints=pd.DataFrame(domain_full.drop(columns='geometry'))
        domain_centpoints_gdf= gpd.GeoDataFrame(domain_centpoints,
            crs='EPSG:4326',
            # in zip() use the coordinate dimentions in form (x coordinate, y coordinate)
            geometry=[shapely.geometry.Point(lonlat) for lonlat in zip(domain_centpoints.LON, domain_centpoints.LAT)])
        domain_centpoints_sub=gpd.clip(domain_centpoints_gdf,state_sub.to_crs(4326))
        
        domain_sub=domain_full.loc[domain_full['ij_key'].isin(domain_centpoints_sub['ij_key']),:]
        
        domain_sub.to_file(state_domain_fn,driver="GPKG")
    
    return state_sub,counties_sub,domain_sub
    
#end make_state_domain

# Read in the drilling, well and flare datasets
def initialize_data(inputs,state,year,month):
    warning=None
    
    # Read in state production totals for monthly scalings
    prod_total_fn = inputs.prod_summary_fn_base.replace('[ST]',state).replace('YYYY',str(year))
    if exists(prod_total_fn):
        state_prod_totals = pd.read_csv(prod_total_fn)
        prod_exist = True
    else: 
        prod_exist = False
    
    
    ####################
    # Drilling data
    ####################
    
    drill_data_fn = inputs.drill_data_fn_base.replace('[ST]',state)
    if exists(drill_data_fn):
        #read all of the drill data in, select only drills occuring in the month, convert to gdf, normalize by number of days in month
        print('Reading in drill data for: '+ state)
        drill_data_df=pd.read_csv(drill_data_fn)
        drill_data_df.rename(columns={'Surface Hole Latitude (WGS84)':'Lat','Surface Hole Longitude (WGS84)':'Lon'},inplace=True)
        drill_data_df['Year']=drill_data_df['Spud Date'].astype(str).str[0:4]
        drill_data_df['Month']=drill_data_df['Spud Date'].astype(str).str[5:7]
        drill_index = np.logical_and(drill_data_df['Year']==str(year),drill_data_df['Month']==str(month).zfill(2))
        drill_data_df=drill_data_df.loc[drill_index,:].reset_index(drop=True)
        drill_data_gdf=gpd.GeoDataFrame(drill_data_df, 
                                             geometry=gpd.points_from_xy(drill_data_df.Lon, drill_data_df.Lat,crs=4326))
        drill_data_gdf['DrCount']=np.ones(drill_data_gdf['Lon'].shape)
        drill_data_gdf['DrFrac']=(np.ones(drill_data_gdf['Lon'].shape))/(drill_data_gdf['DrCount'].sum())
    else:
        print('No drill data found for: '+ state)
        drill_data_gdf=None
    
    ####################
    # Production data
    ####################
    well_data_fn= inputs.well_data_fn_base.replace('[ST]',state).replace('YYYY',str(year)).replace('##',str(month).zfill(2))
    if exists(well_data_fn):
        print('Reading in well production data for: '+ state)
        # well production data read in as per day totals
        well_data_df=pd.read_csv(well_data_fn)
        well_data_gdf=gpd.GeoDataFrame(well_data_df, 
                                             geometry=gpd.points_from_xy(well_data_df.Lon, well_data_df.Lat,crs=4326))
    else:
        print('No well production found for: '+ state)
        well_data_gdf = None
    
    ####################
    # Flare data
    ####################
    
    
    flare_fn = inputs.flare_data_fn_base.replace('YYYY',str(year))
    # if the flare file exists scale based on oil production within the year
    # else scale from the most recent year
    if (exists(flare_fn) and prod_exist):
        print('Reading in flare data')
        flare_data_df=pd.read_csv(flare_fn)
        flare_data_gdf=gpd.GeoDataFrame(flare_data_df, 
                                             geometry=gpd.points_from_xy(flare_data_df.lon, flare_data_df.lat,crs=4326))
        # Generate scaling factor for month and scale flare data
        denom = state_prod_totals.loc[state_prod_totals['Month']==0,'Gprod'].values[0]
        if denom !=0:
            sf = state_prod_totals.loc[state_prod_totals['Month']==month,'Gprod'].values[0]/denom
        else: 
            sf = 1
        flare_data_gdf['flare_NOx']=flare_data_gdf['flare_NOx']*sf # Flare nox in mt/d
        # Calculate CO2 using BTU/yr and EPA natural gas CO2 EF from: https://www.epa.gov/system/files/documents/2022-04/ghg_emission_factors_hub.pdf
        # CO2 EF is 53.06 kg/ mmBtu and convert so flare_CO2 also in mt/d
        #If btu/yr in columns use this, otherwise convert m3/yr to btu/yr using factor from other files (26697.8158506574000 btu/m3) since colby file format not consistent
        if 'btu/yr' in flare_data_gdf.columns:
            flare_data_gdf['flare_CO2']=flare_data_gdf['btu/yr']*sf*53.06/10**6/1000/365
        else:
            flare_data_gdf['flare_CO2']=(flare_data_gdf['m3/yr']*26697.8158506574000)*sf*53.06/10**6/1000/365
    elif prod_exist:
        # list flare filenames and find the one with the closest year 
        flare_filenames_list = [basename(z) for z in glob.glob(inputs.flare_data_fn_base.replace('YYYY','*'))]
        flare_years = [int(z[10:14]) for z in flare_filenames_list]
        offset = [abs(int(z[10:14])-year) for z in flare_filenames_list]
        flare_year_prior = flare_years[offset.index(min(offset))]
        flare_fn_prior = inputs.flare_data_fn_base.replace('YYYY',str(flare_year_prior))
        flare_data_df=pd.read_csv(flare_fn_prior)
        flare_data_gdf=gpd.GeoDataFrame(flare_data_df, 
                                             geometry=gpd.points_from_xy(flare_data_df.lon, flare_data_df.lat,crs=4326))
        # Generate scaling factor for month and scale flare data
        state_prod_totals_prior = pd.read_csv(inputs.prod_summary_fn_base.replace('[ST]',state).replace('YYYY',str(flare_year_prior)))
        denom = state_prod_totals_prior.loc[state_prod_totals_prior['Month']==0,'Gprod'].values[0]
        if denom !=0:
            sf = state_prod_totals.loc[state_prod_totals['Month']==month,'Gprod'].values[0]/denom
        else: 
            sf = 1
        flare_data_gdf['flare_NOx']=flare_data_gdf['flare_NOx']*sf
        # Calculate CO2 using BTU/yr and EPA natural gas CO2 EF from: https://www.epa.gov/system/files/documents/2022-04/ghg_emission_factors_hub.pdf
        # CO2 EF is 53.06 kg/ mmBtu and convert so flare_CO2 also in mt/d
        #If btu/yr in columns use this, otherwise convert m3/yr to btu/yr using factor from other files (26697.8158506574000 btu/m3) since colby file format not consistent
        if 'btu/yr' in flare_data_gdf.columns:
            flare_data_gdf['flare_CO2']=flare_data_gdf['btu/yr']*sf*53.06/10**6/1000/365
        else:
            flare_data_gdf['flare_CO2']=(flare_data_gdf['m3/yr']*26697.8158506574000)*sf*53.06/10**6/1000/365
        
        warning_message = 'WARNING: Flare data for [State: '+ state+ ' Year: ' + str(year) + '] is not present. Using scaled ' +str(flare_year_prior)+' data.'
        print(warning_message)
        warning = {'Code':1,'Message': warning_message}
    else: 
        flare_data_gdf = None
    # End if 
        
        
        
    return drill_data_gdf,well_data_gdf,flare_data_gdf, warning
# end initialize_data

# Joins the drilling, well and flare datasets with the domain
def join_with_domain(drill_data_gdf,well_data_gdf,flare_data_gdf,domain_sub,state):
    domain=domain_sub[['ij_key','Row','Col','geometry']]
    
    print('Joining data with domain for: '+state)
    if (drill_data_gdf is not None):
        drill_sjoin= gpd.sjoin(drill_data_gdf,domain,how='inner').drop(columns='index_right')
        # remove values outside the domain by removing null ij_key
        drill_sjoin_out = drill_sjoin.loc[~drill_sjoin['ij_key'].isnull(),:]
    else:
        drill_sjoin_out = None
        
    if (well_data_gdf is not None):
        well_sjoin= gpd.sjoin(well_data_gdf,domain,how='inner').drop(columns='index_right')
        # remove values outside the domain by removing null ij_key
        well_sjoin_out = well_sjoin.loc[~well_sjoin['ij_key'].isnull(),:]
    else:
        well_sjoin_out = None
        
    if (flare_data_gdf is not None):
        flare_sjoin= gpd.sjoin(flare_data_gdf,domain,how='inner').drop(columns='index_right')
        # remove values outside the domain by removing null ij_key
        flare_sjoin_out = flare_sjoin.loc[~flare_sjoin['ij_key'].isnull(),:]
    else:
        flare_sjoin_out = None
    
    return drill_sjoin_out,well_sjoin_out,flare_sjoin_out
# end join_with_domain
    
# Dissolves the points contained within drill, well and flare data
def dissolve_data(drill_data_gdf,well_data_gdf,flare_data_gdf,state):
    
    # Fill any NA values prior to dissolving
    if (drill_data_gdf is not None):
        drill_dissolve=drill_data_gdf[['ij_key','DrCount','DrFrac','geometry']].copy(deep=True)
        drill_dissolve['DrCount'].fillna(0,inplace=True), drill_dissolve['DrFrac'].fillna(0,inplace=True)
    # end if
    
    if (well_data_gdf is not None):
        
        well_dissolve=well_data_gdf[['ij_key','Oprod','Gprod','DyCO2','LcCO2','geometry']].copy(deep=True)
        well_dissolve['Oprod'].fillna(0,inplace=True)
        well_dissolve['Gprod'].fillna(0,inplace=True)
        well_dissolve['DyCO2'].fillna(0,inplace=True)
        well_dissolve['LcCO2'].fillna(0,inplace=True)
    # end if
    
    if (flare_data_gdf is not None):
        flare_dissolve=flare_data_gdf[['ij_key','flare_NOx','flare_CO2','geometry']].copy(deep=True)
        flare_dissolve['flare_NOx'].fillna(0,inplace=True)
        flare_dissolve['flare_CO2'].fillna(0,inplace=True)
    # end if
    
    # Dissolve datasets
    print('Dissolving data for: '+state)
    if (drill_data_gdf is not None):
        drill_dissolve_out= drill_dissolve.dissolve(by='ij_key',aggfunc='sum').reset_index()
    else:
        drill_dissolve_out=None
    # end if
    if (well_data_gdf is not None):
        well_dissolve_out= well_dissolve.dissolve(by='ij_key',aggfunc='sum').reset_index()
    else:
        well_dissolve_out=None
    # end if
    if (flare_data_gdf is not None):
        flare_dissolve_out= flare_dissolve.dissolve(by='ij_key',aggfunc='sum').reset_index()
    else:
        flare_dissolve_out=None
    # end if
    
    return drill_dissolve_out,well_dissolve_out,flare_dissolve_out
# end join_with_domain

# Join data onto domain by ij_key
# if the dataset does not exist or is empty then fill those spaces with zeros 
def join_onto_domain(drill_dissolve, well_dissolve, flare_dissolve,domain,state):
    print('Joining data onto domain for state: '+state)
    
    # Initialize final data in case the dissolve doesn't exist
    final_data = domain.copy(deep=True)
    
    # If the dataset is none or is empty then fill corresponding columns with zeros 

    if drill_dissolve is not None and not drill_dissolve.empty:
        final_data=domain.join(drill_dissolve.drop(columns='geometry').set_index('ij_key'),on='ij_key')
    else:
        final_data['DrCount']=np.zeros(final_data['ij_key'].shape)
        final_data['DrFrac']=np.zeros(final_data['ij_key'].shape)
    # end if
        

    if well_dissolve is not None and not well_dissolve.empty:
        final_data=final_data.join(well_dissolve.drop(columns='geometry').set_index('ij_key'),on='ij_key')
    else:
        print('setting to zeros')
        final_data['Oprod']=np.zeros(final_data['ij_key'].shape)
        final_data['Gprod']=np.zeros(final_data['ij_key'].shape)
        final_data['DyCO2']=np.zeros(final_data['ij_key'].shape)
        final_data['LcCO2']=np.zeros(final_data['ij_key'].shape)
    # end if
    
    if flare_dissolve is not None and not flare_dissolve.empty:
        final_data=final_data.join(flare_dissolve.drop(columns='geometry').set_index('ij_key'),on='ij_key')
    else:
        final_data['flare_NOx']=np.zeros(final_data['ij_key'].shape)
        final_data['flare_CO2']=np.zeros(final_data['ij_key'].shape)
    # end if
    
    # Fill all NA values with zeros
    final_data['DrCount'].fillna(0,inplace=True)
    final_data['DrFrac'].fillna(0,inplace=True)
    final_data['Oprod'].fillna(0,inplace=True)
    final_data['Gprod'].fillna(0,inplace=True)
    final_data['DyCO2'].fillna(0,inplace=True)
    final_data['LcCO2'].fillna(0,inplace=True)
    final_data['flare_NOx'].fillna(0,inplace=True)
    final_data['flare_CO2'].fillna(0,inplace=True)
    return final_data
# end join onto domain

def add_columns(final_data):
    final_data=final_data.copy(deep=True)
    
    final_data['WcDemand']=np.zeros(final_data['ij_key'].shape)
    final_data['LfElectric']=np.zeros(final_data['ij_key'].shape)
    final_data['LfDemand']=np.zeros(final_data['ij_key'].shape)
    final_data['WetGas']=np.zeros(final_data['ij_key'].shape)
    final_data['WcControll']=np.zeros(final_data['ij_key'].shape)
    final_data['LcControll']=np.zeros(final_data['ij_key'].shape)
    final_data['SUM_HtdryC']=np.zeros(final_data['ij_key'].shape)
    final_data['SUM_HtwetC']=np.zeros(final_data['ij_key'].shape)
    final_data['Sum_WcCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['SUM_LfCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['LeasEngCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['DyFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['LfFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['HtFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['LcFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['WcFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['DrFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['TotFCO2']=np.zeros(final_data['ij_key'].shape)
    final_data['LcFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['WcFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['LfFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['HtFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['DyFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['DrFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['TotFNOx']=np.zeros(final_data['ij_key'].shape)
    final_data['NoppFNOx']=np.zeros(final_data['ij_key'].shape)
    
    final_data.loc[final_data['Oprod']>0,'WetGas'] = 1

    return final_data
    
# end add_columns

def add_spatial(final_data,domain_sub,counties_sub,state):
    # read in spatial_info workbook
    spatial_info=pd.read_excel(inputs.spatial_info_fn,'Spatial_info')
    state_spatial=spatial_info.loc[spatial_info['STATE_ABBR']==state,:]
    
    final_data=final_data.copy(deep=True)
    
    print('Setting spatial info for state: '+state)
    # Calculate centerpoints of domain for selecting which grid cells are in the county
    domain_centpoints=pd.DataFrame(domain_sub.drop(columns='geometry'))
    domain_centpoints_gdf= gpd.GeoDataFrame(domain_centpoints,
        crs='EPSG:4326',
        # in zip() use the coordinate dimentions in form (x coordinate, y coordinate)
        geometry=[shapely.geometry.Point(lonlat) for lonlat in zip(domain_centpoints.LON, domain_centpoints.LAT)])
    
    # Loop over counties that have info in spatial_info and set info in final_data based on this 
    for county in state_spatial['COUNTY_NAME'].unique():
        county_shp=counties_sub.loc[counties_sub['NAME']==county,:]
        county_centpoints=gpd.clip(domain_centpoints_gdf,county_shp.to_crs(4326))
        county_inds=final_data['ij_key'].isin(county_centpoints['ij_key'])
        final_data.loc[county_inds,'LfDemand']=state_spatial.loc[state_spatial['COUNTY_NAME']==county,'LfD'].values[0]
        final_data.loc[county_inds,'LfElectric']=state_spatial.loc[state_spatial['COUNTY_NAME']==county,'LfE'].values[0]
        final_data.loc[county_inds,'WcDemand']=state_spatial.loc[state_spatial['COUNTY_NAME']==county,'WcD'].values[0]
        final_data.loc[county_inds,'WcControll']=state_spatial.loc[state_spatial['COUNTY_NAME']==county,'WcC'].values[0]
        final_data.loc[county_inds,'LcControll']=state_spatial.loc[state_spatial['COUNTY_NAME']==county,'LcC'].values[0]
    
    return final_data
    
# end add_spatial

def add_eia(final_data,state,year,month,scale_eia,drill_counts):
    final_data=final_data.copy(deep=True)
    
    # read in EIA data
    eia_data=pd.read_excel(inputs.eia_info_fn,'eia')
    eia_data_drilling = eia_data.dropna(subset=['Diesel'])
    eia_data_lease = eia_data.dropna(subset=['Lease CO2'])
    state_eia_drilling=eia_data_drilling.loc[np.logical_and(eia_data_drilling['State']==state, eia_data_drilling['Year']==year),:]
    state_eia_lease=eia_data_lease.loc[np.logical_and(eia_data_lease['State']==state, eia_data_lease['Year']==year),:]
    
    # Read in state production totals for monthly scalings
    prod_total_fn = inputs.prod_summary_fn_base.replace('[ST]',state).replace('YYYY',str(year))
    if exists(prod_total_fn):
        state_prod_totals = pd.read_csv(prod_total_fn)
        prod_exist = True
    else: 
        prod_exist = False
    
    warning=None
    no_eia_drilling=False
    if state_eia_drilling.empty:
        i=0
        while state_eia_drilling.empty:
            i=i+1
            state_eia_drilling=eia_data_drilling.loc[np.logical_and(eia_data_drilling['State']==state, eia_data_drilling['Year']==year-i),:]
            if i > 50:
                no_eia_drilling = True# if too many searches without anything then set values such that eia_lease and offdsl are zero
                break
        # End while
        
        # Calculate drilling scaling factor 
        if drill_counts is not None:
            num_days = calendar.monthrange(year,month)[1]
            drill_count_annprior = drill_counts.loc[drill_counts['Year']==str(year-i),'DrCount'] 
            drill_count_month = drill_counts.loc[np.logical_and(drill_counts['Year']==str(year),drill_counts['Month']==str(month).zfill(2)),'DrCount']
            
            if drill_count_annprior.empty or drill_count_month.empty: 
                drill_sf=0
                warning_message = 'WARNING: Drilling Info for [State: '+ state+ ' Year: ' + str(year) + ' or ' + ' Year: '+ str(year-i) + ' Month: '+str(month)+'] is not present.  Setting drill diesel scaling to zero.'
                print(warning_message)
                if warning:
                    warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
                else:
                    warning = {'Code':1,'Message':  warning_message}
            else:
                drill_count_annprior = drill_count_annprior.sum()/365# annual average drill count per day from prior year
                drill_count_month = drill_count_month.values[0]/num_days # monthly average drill count per day from current month
                drill_sf = drill_count_month/drill_count_annprior
                
                warning_message = 'WARNING: EIA Info for [State: '+ state+ ' Year: ' + str(year) + '] is not present. Using ' +str(year-i) +' data for drilling diesel, scaling based on drill count.'
                print(warning_message)
                if warning:
                    warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
                else:
                    warning = {'Code':1,'Message':  warning_message}
        else:
            drill_sf = 0
            warning_message = 'WARNING: drill file for State: '+ state +' is not present. Setting drill diesel scaling to zero.'
            print(warning_message)
            if warning:
                warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
            else:
                warning = {'Code':1,'Message':  warning_message}
    else:
        # Calculate drilling scaling factor 
        if drill_counts is not None:
            num_days = calendar.monthrange(year,month)[1]
            drill_count_ann = drill_counts.loc[drill_counts['Year']==str(year),'DrCount'] 
            drill_count_month = drill_counts.loc[np.logical_and(drill_counts['Year']==str(year),drill_counts['Month']==str(month).zfill(2)),'DrCount']
            
            if drill_count_ann.empty or drill_count_month.empty: 
                drill_sf=0
                warning_message = 'WARNING: Drilling Info for [State: '+ state+ ' Year: ' + str(year) + ' Month: '+str(month)+'] is not present.  Setting drill diesel scaling to zero.'
                print(warning_message)
                if warning:
                    warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
                else:
                    warning = {'Code':1,'Message':  warning_message}
            else:
                drill_count_ann = drill_count_ann.sum()/365# annual average drill count per day from prior year
                drill_count_month = drill_count_month.values[0]/num_days # monthly average drill count per day from current month
                drill_sf = drill_count_month/drill_count_ann

        else:
            drill_sf = 0
            warning_message = 'WARNING: drill file for State: '+ state +' is not present. Setting drill diesel scaling to zero.'
            print(warning_message)
            if warning:
                warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
            else:
                warning = {'Code':1,'Message':  warning_message}
        #end if
    
    if no_eia_drilling:
        offdsl = 0 
    else:
        offdsl=state_eia_drilling['Diesel'].values[0]/365*drill_sf
    # end if
    
    # If the EIA info for the year doesn't exist then use the previous years data but return a warning 
    # Use the gas production data to scale the eia lease gas
    no_eia_lease=False
    if state_eia_lease.empty:
        i=0
        while state_eia_lease.empty:
            i=i+1
            state_eia_lease=eia_data_lease.loc[np.logical_and(eia_data_lease['State']==state, eia_data_lease['Year']==year-i),:]
            if i > 50:
                no_eia_lease = True# if too many searches without anything then set values such that eia_lease and offdsl are zero
                break
        # End while
        
        if (~no_eia_lease and prod_exist):
            state_prod_totals_prior = pd.read_csv(inputs.prod_summary_fn_base.replace('[ST]',state).replace('YYYY',str(year-i)))
            denom = state_prod_totals_prior.loc[state_prod_totals_prior['Month']==0,'Gprod'].values[0]
            if denom !=0:
                prod_sf = (state_prod_totals.loc[state_prod_totals['Month']==month,'Gprod'].values[0]/
                           denom)
            else:
                prod_sf = 1
            warning_message = 'WARNING: EIA Info for [State: '+ state+ ' Year: ' + str(year) + '] is not present. Using ' +str(year-i) +' data for lease gas, scaling based on gas production.'
            print(warning_message)
            if warning:
                warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
            else:
                warning = {'Code':1,'Message':  warning_message}
        elif ~prod_exist: 
            prod_sf = 0
            warning_message = 'WARNING: Production for [State: '+ state+ ' Year: ' + str(year) + '] is not present. Setting lease gas scale factor to zero.'
            print(warning_message)
            if warning:
                warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
            else:
                warning = {'Code':1,'Message':  warning_message}

    elif prod_exist:
        warning = None
        denom =state_prod_totals.loc[state_prod_totals['Month']==0,'Gprod'].values[0]
        if denom !=0:
            prod_sf = state_prod_totals.loc[state_prod_totals['Month']==month,'Gprod'].values[0]/denom
        else:
            prod_sf = 1
        no_eia_lease = False
        
    else: 
        prod_sf = 0
        warning_message = 'WARNING: Production for [State: '+ state+ ' Year: ' + str(year) + '] is not present. Setting lease gas factor to zero.'
        print(warning_message)
        if warning:
            warning = {'Code':1,'Message': warning['Message']+ ' \n' + warning_message}
        else:
            warning = {'Code':1,'Message':  warning_message}
    
    # Print the scaling factors used
    message = 'For [State: '+ state+ ' Year: ' + str(year) + ' Month: ' +str(month)+ '] Drill scaling = ' + str(drill_sf) +' Lease Gas scaling = '+str(prod_sf)
    print(message)
    
    if no_eia_lease:
        eia_lease = 0
    else:
        eia_lease=state_eia_lease['Lease CO2'].values[0]/365*prod_sf
    # end if
    
    
    # Begin Calculating fields
    print('Calculating new fields for state: '+state)
    final_data['SUM_HtdryC']=final_data['Gprod']*0.702159198152654*(1-final_data['WetGas'])
    final_data['SUM_HtwetC']=final_data['Gprod']*final_data['WetGas']*1.99634880053995
    final_data['Sum_WcCO2']=final_data['Gprod']*2.61048580797312*final_data['WcDemand']
    final_data['SUM_LfCO2']=final_data['Oprod']*8.29883410616054*final_data['LfDemand']*(1-final_data['LfElectric'])
    final_data['LeasEngCO2']=final_data['DyCO2']+final_data['LcCO2']+final_data['SUM_HtdryC']+final_data['SUM_HtwetC']+final_data['Sum_WcCO2']+final_data['SUM_LfCO2']
    
    LeaseEngCO2_sum=final_data['LeasEngCO2'].sum() # Units of per day
    if LeaseEngCO2_sum==0:
        LeaseEngCO2_sum=1
        
    # CO2 fields
    final_data['LfFCO2']=final_data['SUM_LfCO2']/LeaseEngCO2_sum*eia_lease/1000 # divide by 1000 so in mt/d
    final_data['DyFCO2']=final_data['DyCO2']/LeaseEngCO2_sum*eia_lease/1000 # divide by 1000 so in mt/d
    final_data['LcFCO2']=final_data['LcCO2']/LeaseEngCO2_sum*eia_lease/1000 # divide by 1000 so in mt/d
    final_data['WcFCO2']=final_data['Sum_WcCO2']/LeaseEngCO2_sum*eia_lease/1000 # divide by 1000 so in mt/d
    final_data['HtFCO2']=(final_data['SUM_HtdryC']+final_data['SUM_HtwetC'])/LeaseEngCO2_sum*eia_lease/1000 # divide by 1000 so in mt/d
    final_data['DrFCO2']=final_data['DrFrac']*10.16*1000*offdsl/1000 # offdsl multiplied by 1000 because in 1000's of gallons need to convert to gallons, divide by 1000 so in mt/d
    final_data['TotFCO2']=(final_data['DyFCO2']+final_data['LfFCO2']+final_data['HtFCO2']+final_data['LcFCO2']+final_data['WcFCO2']+final_data['DrFCO2']) +final_data['flare_CO2'] # This in mt/day 
    
    # NOx fields, This all in mt/d 
    final_data['DrFNOx']=final_data['DrFCO2']*11.5/1000 # divide by 1000 so ratio is in kg/kg and final number in mt/d
    final_data['DyFNOx']=final_data['DyFCO2']*0.76/1000 # divide by 1000 so ratio is in kg/kg and final number in mt/d
    final_data['HtFNOx']=final_data['HtFCO2']*0.83/1000 # divide by 1000 so ratio is in kg/kg and final number in mt/d
    final_data['LfFNOx']=final_data['LfFCO2']*20.5/1000 # divide by 1000 so ratio is in kg/kg and final number in mt/d
    final_data['LcFNOx']=(final_data['LcFCO2']*1.57*final_data['LcControll'] + (1-final_data['LcControll'])*4.24*final_data['LcFCO2'])/1000 # divide by 1000 so ratio is in kg/kg and final number in mt/d
    final_data['WcFNOx']=(final_data['WcFCO2']*35.35*(1-final_data['WcControll']) + final_data['WcControll']*4.33*final_data['WcFCO2'])/1000 # divide by 1000 so ratio is in kg/kg and final number in mt/d
    final_data['NoppFNOx']=(final_data['WcFNOx']+final_data['LcFNOx']+final_data['DyFNOx']+final_data['DrFNOx']+final_data['HtFNOx']+final_data['LfFNOx']) +final_data['flare_NOx']
    
    return final_data, warning
# end add_eia

# Add all VOC species based on basin specific ratios
def add_voc_species(final_data,domain_sub,basins_shp,play_shp,state,year,month,state_ratios):
    
    final_data=final_data.copy(deep=True)
    
    print('Add VOC species to state: '+state)
    
    # Add all species as zero values
    final_data['CH4']=np.zeros(final_data['ij_key'].shape)
    final_data['ALK1']=np.zeros(final_data['ij_key'].shape)
    final_data['ALK4']=np.zeros(final_data['ij_key'].shape)
    final_data['ALK5']=np.zeros(final_data['ij_key'].shape)
    final_data['ARO1']=np.zeros(final_data['ij_key'].shape)
    final_data['BENZENE']=np.zeros(final_data['ij_key'].shape)
    final_data['BUTANES']=np.zeros(final_data['ij_key'].shape)
    final_data['PENTANES']=np.zeros(final_data['ij_key'].shape)
    final_data['TOLUENE']=np.zeros(final_data['ij_key'].shape)
    final_data['M_XYLENE']=np.zeros(final_data['ij_key'].shape)
    final_data['O_XYLENE']=np.zeros(final_data['ij_key'].shape)
    final_data['P_XYLENE']=np.zeros(final_data['ij_key'].shape)
    final_data['PROPANE']=np.zeros(final_data['ij_key'].shape)
    final_data['NMHC']=np.zeros(final_data['ij_key'].shape)
    
    final_data['CH4_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['ALK1_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['ALK4_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['ALK5_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['ARO1_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['BENZENE_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['BUTANES_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PENTANES_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['TOLUENE_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['M_XYLENE_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['O_XYLENE_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['P_XYLENE_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PROPANE_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['NMHC_Dr']=np.zeros(final_data['ij_key'].shape)
    
    # Slopes of correlations with methane
    m = {'SNMHC':1.26,
    'SBenzene':0.0030693898/78.11,
    'SBenzene_1ethyl':0.000450778/106.165,
    'SButane_iso':0.08/58.12,
    'SButane_n':0.27/58.12,
    'SCyclohexane':0.0045/84.16,
    'SCyclohexane_1methyl':0.0045/98.186,
    'SEthane':0.34/30.07,
    #SEthene=str()
    #SEthyne=str()
    'SHeptane_n':0.008/100.21,
    'SHexane_n':0.024/86.18,
    'SOctane_n':0.002/114.23,
    'SPentane_iso':0.062/72.15,
    'SPentane_n':0.092/72.15,
    'SPropane':0.37/44.1,
    'SToluene':0.003383352/92.14,
    'SMP_XYLENE':0.001171747/106.16,
    'SO_XYLENE':0.000372077/106.16}
    
    # Read in basin ratios workbook
    basin_fn = (inputs.basin_info_fn_base)
    
    # If basin ratios do not exist for this year, search for other years

    basin_data_ofrac=pd.read_excel(basin_fn,'Ofrac')
    basin_data_ch4_ratio=pd.read_excel(basin_fn,'CH4_ratio')
    
    basin_names = basin_data_ch4_ratio['Basin name'].unique()
    
    # read in the names of the plays for methane trends
    play_names = play_shp['Play_Name'].unique()
    
    # Read in play methane scaling data
    play_trends = pd.read_excel(inputs.play_methane_trends_fn,'Sheet1')
    
    # Find US methane trend scaling
    US_methane_trends = play_trends.loc[play_trends['Play_Name']=='US',:]
    US_trends_year_max = US_methane_trends['Year'].max()
    US_trends_year_min = US_methane_trends['Year'].min()
    if year > US_trends_year_max:
        US_intensity_ratio = US_methane_trends.loc[US_methane_trends['Year']==US_trends_year_max,'scaling_rel_to_2015'].values[0]
    elif year < US_trends_year_min:
        US_intensity_ratio = US_methane_trends.loc[US_methane_trends['Year']==US_trends_year_min,'scaling_rel_to_2015'].values[0]
    else:
        US_intensity_ratio = US_methane_trends.loc[US_methane_trends['Year']==year,'scaling_rel_to_2015'].values[0]
    
    domain_centpoints=pd.DataFrame(domain_sub.drop(columns='geometry'))
    domain_centpoints_gdf= gpd.GeoDataFrame(domain_centpoints,
        crs='EPSG:4326',
        # in zip() use the coordinate dimentions in form (x coordinate, y coordinate)
        geometry=[shapely.geometry.Point(lonlat) for lonlat in zip(domain_centpoints.LON, domain_centpoints.LAT)])
    
    # Loop over basins for Ofrac and CH4 ratio
    for basin in basin_names:
        basin_shp_sub=basins_shp.loc[basins_shp['BASIN_NAME']==basin,:]
        basin_centpoints=gpd.clip(domain_centpoints_gdf,basin_shp_sub.to_crs(4326))
        basin_inds=final_data['ij_key'].isin(basin_centpoints['ij_key'])
        
        ch4_ratio=basin_data_ch4_ratio.loc[basin_data_ch4_ratio['Basin name']==basin,'CH4_ratio'].values[0]
        Ofrac=basin_data_ofrac.loc[np.logical_and(basin_data_ofrac['Year']==year,basin_data_ofrac['Month']==month),basin].values[0]
        
        # First do for all with US methane scaling scaling
        if any(basin_inds):
            print('For Basin: ' + basin + ' US Methane Intensity Trend Scaling: ' + str(US_intensity_ratio) )
            final_data = calc_ch4(final_data,basin_inds,ch4_ratio,state,state_ratios,US_intensity_ratio)
        else: 
            continue
        
        # Loop over different EIA plays to apply methane trends relative to 2015
        for play in play_names:
            play_shp_sub=play_shp.loc[play_shp['Play_Name']==play,:]
            play_centpoints=gpd.clip(domain_centpoints_gdf,play_shp_sub.to_crs(4326))
            #if basin == 'Denver':
            #    breakpoint()
            basin_play_inds=final_data['ij_key'].isin(basin_centpoints['ij_key']) & final_data['ij_key'].isin(play_centpoints['ij_key'])
            

            # if any inds are in the combo of basin and play then recalculate with the appropriate play specific trend, otherwise cycle
            if any(basin_play_inds):
                
                # Find play methane trend scaling
                play_methane_trends = play_trends.loc[play_trends['Play_Name']==play,:]
                play_methane_trends_year_max = play_methane_trends['Year'].max()
                play_methane_trends_year_min = play_methane_trends['Year'].min()
                if year > play_methane_trends_year_max:
                    play_intensity_ratio = play_methane_trends.loc[play_methane_trends['Year']==play_methane_trends_year_max,'scaling_rel_to_2015'].values[0]
                elif year < play_methane_trends_year_min:
                    play_intensity_ratio = play_methane_trends.loc[play_methane_trends['Year']==play_methane_trends_year_min,'scaling_rel_to_2015'].values[0]
                else:
                    play_intensity_ratio = play_methane_trends.loc[play_methane_trends['Year']==year,'scaling_rel_to_2015'].values[0]
                
                print('Found Play specific methane trend within basin, applying now: ' )
                print('For Basin: ' + basin + ' Play: ' + play + ' Play Methane Intensity Trend Scaling: ' + str(play_intensity_ratio) )
                final_data = calc_ch4(final_data,basin_play_inds,ch4_ratio,state,state_ratios,play_intensity_ratio)
            else:
                continue
        
        
        # if any inds are in the state domain then proceed, otherwise cycle
        if any(basin_inds):
            print('Calculating VOC species other than methane: ' )
            print('For Basin: ' + basin + ' Oil Fraction : ' +str(Ofrac))
            final_data = calc_voc_species(final_data,m,basin_inds,Ofrac)
        else:
            continue
    # end for loop over basins
    
    warning=None
    return final_data, warning
# end add_VOC_species

# Calculates CH4 for each basin, adding methane trends for play boundaries
def calc_ch4(final_data,inds,ch4_ratio,state,state_ratios,intensity_ratio):
    final_data=final_data.copy(deep=True)
    # calculate all VOC species
    
    ProdNOx = final_data['NoppFNOx'].sum() -final_data['DrFNOx'].sum()# total NOx this year
    TotNOx = final_data['NoppFNOx'].sum()# total NOx this year
    ratio_2015 =state_ratios.loc[state_ratios['State']==state,'Ratio']
    if ratio_2015.empty:
        ratio_2015 = 0 # set this to zero so that ratio gets set to 1 down lower
    else:
        ratio_2015 = state_ratios.loc[state_ratios['State']==state,'Ratio'].values[0]
    
    # add if statements to avoid NaN values
    if TotNOx == 0:
        ratio=1
    else:
        ratio = ratio_2015*(ProdNOx/TotNOx)
        if ratio==0:
            ratio=1

    
    print('NOx scaling Ratio: '+ str(ratio))
    
    # final_data['CH4']=final_data['NoppFNOx']*1000000*1024/677*750*1168/46/ch4_ratio # Ratio here is Total NOx 2015/ProductionNOx 2015 * Production NOx this year/ Total NOx this year
    # CH4 in mol/day
    # All other VOCs below also in mol/day because S${VOC} includes the molecular weight of each species in denominator 
    final_data.loc[inds,'CH4']=final_data.loc[inds,'NoppFNOx']*1000000*ratio/46/ch4_ratio*intensity_ratio # Ratio here is Total NOx 2015/ProductionNOx 2015 * Production NOx this year/ Total NOx this year
    
    final_data.loc[inds,'CH4_Dr']=(final_data.loc[inds,'DrFNOx'])*1000000*ratio/46/ch4_ratio*intensity_ratio # Ratio here is Total NOx 2015/ProductionNOx 2015 * Production NOx this year/ Total NOx this year
    
    return final_data
    
# calculates VOC species for each basin
def calc_voc_species(final_data,m,basin_inds,Ofrac):
    
    final_data=final_data.copy(deep=True)
    # calculate all VOC species
    
    final_data.loc[basin_inds,'NMHC']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SNMHC'] 
    final_data.loc[basin_inds,'BENZENE']=final_data.loc[basin_inds,'CH4']*16*m['SBenzene']
    final_data.loc[basin_inds,'ARO1']=final_data.loc[basin_inds,'CH4']*16*m['SBenzene_1ethyl']
    final_data.loc[basin_inds,'BUTANES']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SButane_iso']+final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SButane_n']
    final_data.loc[basin_inds,'ALK1']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SEthane']
    final_data.loc[basin_inds,'ALK4']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SHeptane_n']+final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SHexane_n']
    final_data.loc[basin_inds,'ALK5']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SCyclohexane']+final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SCyclohexane_1methyl'] + final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SOctane_n']
    final_data.loc[basin_inds,'PENTANES']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SPentane_iso']+final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SPentane_n']
    final_data.loc[basin_inds,'PROPANE']=final_data.loc[basin_inds,'CH4']*16*Ofrac*m['SPropane']
    final_data.loc[basin_inds,'TOLUENE']=final_data.loc[basin_inds,'CH4']*16*m['SToluene']
    final_data.loc[basin_inds,'M_XYLENE']=final_data.loc[basin_inds,'CH4']*16*0.5*m['SMP_XYLENE']
    final_data.loc[basin_inds,'O_XYLENE']=final_data.loc[basin_inds,'CH4']*16*m['SO_XYLENE']
    final_data.loc[basin_inds,'P_XYLENE']=final_data.loc[basin_inds,'CH4']*16*0.5*m['SMP_XYLENE']
    
    final_data.loc[basin_inds,'NMHC_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SNMHC'] 
    final_data.loc[basin_inds,'BENZENE_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*m['SBenzene']
    final_data.loc[basin_inds,'ARO1_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*m['SBenzene_1ethyl']
    final_data.loc[basin_inds,'BUTANES_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SButane_iso']+final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SButane_n']
    final_data.loc[basin_inds,'ALK1_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SEthane']
    final_data.loc[basin_inds,'ALK4_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SHeptane_n']+final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SHexane_n']
    final_data.loc[basin_inds,'ALK5_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SCyclohexane']+final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SCyclohexane_1methyl'] + final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SOctane_n']
    final_data.loc[basin_inds,'PENTANES_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SPentane_iso']+final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SPentane_n']
    final_data.loc[basin_inds,'PROPANE_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*Ofrac*m['SPropane']
    final_data.loc[basin_inds,'TOLUENE_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*m['SToluene']
    final_data.loc[basin_inds,'M_XYLENE_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*0.5*m['SMP_XYLENE']
    final_data.loc[basin_inds,'O_XYLENE_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*m['SO_XYLENE']
    final_data.loc[basin_inds,'P_XYLENE_Dr']=final_data.loc[basin_inds,'CH4_Dr']*16*0.5*m['SMP_XYLENE']
    
    
    return final_data
# end calc_voc_species

# Add species based on ratio to CO2
def add_CO2_ratio_species(final_data,state,year,month):
    
    final_data=final_data.copy(deep=True)
    
    print('Add species based on ratio to CO2 species to state: '+state)
    
    # Add all species as zero values
    final_data['CO']=np.zeros(final_data['ij_key'].shape)
    final_data['NH3']=np.zeros(final_data['ij_key'].shape)
    final_data['PM10-PRI']=np.zeros(final_data['ij_key'].shape)
    final_data['PM25-PRI']=np.zeros(final_data['ij_key'].shape)
    final_data['SO2']=np.zeros(final_data['ij_key'].shape)
    final_data['ETHENE']=np.zeros(final_data['ij_key'].shape) # also HC07
    final_data['HCHO']=np.zeros(final_data['ij_key'].shape)# also HC14
    final_data['CCHO']=np.zeros(final_data['ij_key'].shape)# also HC15
    final_data['ACETYLENE']=np.zeros(final_data['ij_key'].shape)# also HC37
    final_data['PM01']=np.zeros(final_data['ij_key'].shape)
    final_data['PM02']=np.zeros(final_data['ij_key'].shape)
    final_data['PM03']=np.zeros(final_data['ij_key'].shape)
    final_data['PM04']=np.zeros(final_data['ij_key'].shape)
    final_data['PM05']=np.zeros(final_data['ij_key'].shape)
    final_data['PM06']=np.zeros(final_data['ij_key'].shape)
    final_data['PM07']=np.zeros(final_data['ij_key'].shape)
    final_data['PM08']=np.zeros(final_data['ij_key'].shape)
    final_data['PM09']=np.zeros(final_data['ij_key'].shape)
    final_data['PM10']=np.zeros(final_data['ij_key'].shape)
    final_data['PM11']=np.zeros(final_data['ij_key'].shape)
    final_data['PM12']=np.zeros(final_data['ij_key'].shape)
    final_data['PM13']=np.zeros(final_data['ij_key'].shape)
    final_data['PM14']=np.zeros(final_data['ij_key'].shape)
    final_data['PM15']=np.zeros(final_data['ij_key'].shape)
    final_data['PM16']=np.zeros(final_data['ij_key'].shape)
    final_data['PM17']=np.zeros(final_data['ij_key'].shape)
    final_data['PM18']=np.zeros(final_data['ij_key'].shape)
    final_data['PM19']=np.zeros(final_data['ij_key'].shape)
    
    # Add all species as zero values
    final_data['CO_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['NH3_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM10-PRI_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM25-PRI_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['SO2_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['ETHENE_Dr']=np.zeros(final_data['ij_key'].shape) # also HC07
    final_data['HCHO_Dr']=np.zeros(final_data['ij_key'].shape)# also HC14
    final_data['CCHO_Dr']=np.zeros(final_data['ij_key'].shape)# also HC15
    final_data['ACETYLENE_Dr']=np.zeros(final_data['ij_key'].shape)# also HC37
    final_data['PM01_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM02_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM03_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM04_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM05_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM06_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM07_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM08_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM09_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM10_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM11_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM12_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM13_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM14_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM15_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM16_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM17_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM18_Dr']=np.zeros(final_data['ij_key'].shape)
    final_data['PM19_Dr']=np.zeros(final_data['ij_key'].shape)
    
    # Ratio to CO2 (units differ)
    m = {'CO':0.007763225, # mt/mt CO2
        'NH3':8.79219E-08, # mt/mt CO2
        'PM10-PRI':0.000134296, # mt/mt CO2
        'PM25-PRI':0.000132737, # mt/mt CO2
        'SO2':0.000443601, # mt/mt CO2
        'ETHENE':46.75519631, # moles/mt CO2
        'HCHO':64.17782826, # moles/mt CO2
        'CCHO':28.74263855, # moles/mt CO2
        'ACETYLENE':99.66639729, # moles/mt CO2
        'PM01':4.59285E-05, # mt/mt CO2
        'PM02':1.88694E-05, # mt/mt CO2
        'PM03':2.30749E-06, # mt/mt CO2
        'PM04':4.25218E-05, # mt/mt CO2
        'PM05':6.10859E-06, # mt/mt CO2
        'PM06':1.70135E-05, # mt/mt CO2
        'PM07':4.08147E-06, # mt/mt CO2
        'PM08':7.14883E-07, # mt/mt CO2
        'PM09':9.64527E-07, # mt/mt CO2
        'PM10':1.42892E-06, # mt/mt CO2
        'PM11':2.44962E-06, # mt/mt CO2
        'PM12':0, # mt/mt CO2
        'PM13':3.11933E-05, # mt/mt CO2
        'PM14':3.36417E-07, # mt/mt CO2
        'PM15':0, # mt/mt CO2
        'PM16':2.52743E-06, # mt/mt CO2
        'PM17':0, # mt/mt CO2
        'PM18':0, # mt/mt CO2
        'PM19':2.23268E-06, # mt/mt CO2
    }
    
    # Loop over and set values for total emissions
    # Values for criteria and PM species after this are metric tons/day 
    # Values for VOCs after this are moles/day due to use of m 
    for key in m.keys():
        final_data.loc[:,key]=final_data.loc[:,'TotFCO2']*m[key]
    
    # Loop over and set values for drilling emissions
    # Values for criteria and PM species after this are metric tons/day 
    # Values for VOCs after this are moles/day
    for key in m.keys():
        final_data.loc[:,key+'_Dr']=final_data.loc[:,'DrFCO2']*m[key]      
    
    
    warning=None
    return final_data, warning
# end add_CO2_ratio_species
    
def save_output(inputs, final_data,state,year,month):
    
    if inputs.save_gpkg:
        out_fn = inputs.output_dir + inputs.out_fn_base.replace('[ST]',state).replace('YYYY',str(year)).replace('##',str(month).zfill(2)) + 'gpkg'
        print('Writing output gpkg for state: '+state)
        final_data.to_file(out_fn,driver="GPKG")
    if inputs.save_pkl:
        out_fn = inputs.output_dir + inputs.out_fn_base.replace('[ST]',state).replace('YYYY',str(year)).replace('##',str(month).zfill(2))+ 'pkl.gz'
        print('Writing output pickle for state: '+state)
        compression_opts={'method': 'gzip', 'compresslevel': 1}
        final_data.drop(columns='geometry').to_pickle(out_fn,compression=compression_opts)
        
# end save_output

def write_warning_text(warning_eia, warning_flare, warning_voc, warning):
    # Add warning messages to log
    if (warning_eia is not None):
        if (warning is None):
            warning = {'Code':1,
                        'Message': warning_eia['Message'] + '\n'}
        else:
            warning['Message'] = warning['Message'] + warning_eia['Message'] + '\n'
    elif (warning_flare is not None):
        if (warning is None):
            warning = {'Code':1,
                        'Message': warning_flare['Message'] + '\n'}
        else:
            warning['Message'] = warning['Message'] + warning_flare['Message'] + '\n'
    elif (warning_voc is not None):
        if (warning is None):
            warning = {'Code':1,
                        'Message': warning_voc['Message'] + '\n'}
        else:
            warning['Message'] = warning['Message'] + warning_voc['Message'] + '\n'
    # end if
    return warning
# end write_warning_text

# Calculates the number of wells drilled for each year and each month, returns as df
def calc_drill_counts(inputs,state):
    print('Calculating drilling summary statistics for: '+ state)
    drill_data_fn = inputs.drill_data_fn_base.replace('[ST]',state)
    if exists(drill_data_fn):
        drill_data_df=pd.read_csv(drill_data_fn)
        drill_data_df['Year']=drill_data_df['Spud Date'].astype(str).str[0:4]
        drill_data_df['Month']=drill_data_df['Spud Date'].astype(str).str[5:7]
        drill_data_df['DrCount']=np.ones(drill_data_df['Year'].shape)
        
        drill_summary = drill_data_df[['Year','Month','DrCount']].groupby(by=['Year','Month']).sum().reset_index()
        
    else:
        print('No drill data found for: '+ state)
        drill_summary=None
        
    return drill_summary
# end calc_drill_counts

def main():
    start_time = time.time()
    print('START OF PROGRAM: ')
    print('--------------------------------------------------')

    print('Reading in national shapefiles....')    
    # Read in shapefiles for whole country so it only happens once
    states_full=gpd.read_file(inputs.states_fn)
    counties_full=gpd.read_file(inputs.counties_fn)
    domain_full=gpd.read_file(inputs.domain_fn)
    basins_shp=gpd.read_file(inputs.basins_fn)
    play_shp=gpd.read_file(inputs.plays_fn)
    
    # Read in state ratios
    state_ratios = pd.read_csv(inputs.state_ratios_fn)
    
    warning = None
    
    # Loop over states in list
    for state in inputs.state_list:
        print('Making State Domain files for: ' +state)
        state_sub,counties_sub,domain_sub = make_state_domain(inputs,state,domain_full,counties_full,states_full)
        
        drill_counts = calc_drill_counts(inputs,state)
        # Loop over years in list
        for year in inputs.years:
            # Loop over months in list
            for month in inputs.months: 
                start_state_time = time.time()
                
                print('Working on Emissions for: ' + str(year) + ', Month : ' + str(month) + ', State: ' + state )
                
                
                drill_data_gdf,well_data_gdf,flare_data_gdf, warning_flare = initialize_data(inputs,state,year,month)
                drill_spatial,well_spatial,flare_spatial = join_with_domain(drill_data_gdf,well_data_gdf,flare_data_gdf,domain_sub,state)
                drill_data_dissolve,well_data_dissolve,flare_data_dissolve = dissolve_data(drill_spatial,well_spatial,flare_spatial,state)
                
                domain_out=domain_sub.copy(deep=True)
                final_data = join_onto_domain(drill_data_dissolve, well_data_dissolve, flare_data_dissolve,domain_out,state)
                final_data = add_columns(final_data)
                final_data = add_spatial(final_data,domain_sub,counties_sub,state)
                
                scale_eia = ~(well_data_dissolve is None) # add a flag so that production scaling is 1 if no production data exists
                final_data, warning_eia = add_eia(final_data,state,year,month,scale_eia,drill_counts)
                final_data, warning_voc = add_voc_species(final_data,domain_sub,basins_shp,play_shp,state,year,month,state_ratios)
                final_data, warning_voc = add_CO2_ratio_species(final_data,state,year,month)
                
                warning = write_warning_text(warning_eia, warning_flare, warning_voc, warning)
                
                # print(final_data)
                
                save_output(inputs, final_data,state,year,month)
                
                end_state_time = time.time()
                # print('Time to make files for Year : ' + str(year) + ',  Month: ' + str(m) + ' is ' + str( round(end_state_time-start_state_time ) )  + ' s')
                print('Time elapsed to make files for Year : ' + str(year) + ', Month : ' + str(month) + ', State: ' + state + ' is ' + str( round(end_state_time-start_state_time ) )  + ' s')
            # end Loop over months
        # end Loop over years
    # end for loop over state
    
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
    
    # print('NOTES: ')

    
    

if __name__ == "__main__":
    main()
