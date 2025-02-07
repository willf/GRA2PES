import pandas as pd
import numpy as np
import requests
import calendar

class inputs():
    # https://medium.com/analytics-vidhya/use-python-to-pull-energy-data-from-the-us-department-of-energys-api-11d6f724927e
    
    # Economic indicators API webpage
    # https://www.census.gov/data/developers/data-sets/economic-indicators.html
    # https://www2.census.gov/data/api-documentation/EITS_API_User_Guide_Dec2020.pdf
    
    # Wholesale trade API address
    #https://api.census.gov/data/timeseries/eits/mwts.html
    
    #XML where you can find variable info links
    
    
    key=''
    # API_base_address = 'https://api.eia.gov/series/?api_key='+key+'&series_id='
    API_base_address = ('https://api.census.gov/data/timeseries/eits/mwts'+
                        '?get=time_slot_id,cell_value&data_type_code=DATATYPECODE&category_code=CATCODE&error_data=no&seasonally_adj=no&for=us:*&time=from+2000&key=' + key)
    
    
    # Data codes to loop over
    dc_sales = 'SM'
    dc_inventory = 'IM'
    
    series_id = {'TotalMerchant':'42', # Total MerchantWholesalers,ExceptManufacturers'Sales Branches and Offices
                 'MotorVehicle':'4231',# MotorVehicle &MotorVehicleParts &Supplies
                  'Lumber':'4233', # Lumber &OtherConstructionMaterials
                  'ProfessionalSupplies':'4234',#Professional& CommercialEquipment &Supplies
                  'Metals':'4235', # Metals &Minerals,ExceptPetroleum
                  'Electronics':'4236', #HouseholdAppliances& Electrical& ElectronicGoods
                  'Machinery':'4238', # Machinery,Equipment,& Supplies
                  'Drug':'4242', # Drugs &Druggists'Sundries
                  'Apparel':'4243', # Apparel,PieceGoods, & Notions
                  'Grocery':'4244', # Grocery &RelatedProducts
                  'Chemical':'4246', # Chemicals& AlliedProducts
                  'Petroleum':'4247', # Petroleum& PetroleumProducts
                  'Misc':'4249'
                 } 
    
    ppi_id = {'TotalMerchant':'PCUAWHLTRAWHLTR', # Total MerchantWholesalers,ExceptManufacturers'Sales Branches and Offices
                 'MotorVehicle':'PCU4231--4231--',# MotorVehicle &MotorVehicleParts &Supplies
                 'Lumber':'PCU4233--4233--', # Lumber &OtherConstructionMaterials
                 'ProfessionalSupplies':'PCU4234--4234--',#Professional& CommercialEquipment &Supplies
                 'Metals':'PCU4235--4235--', # Metals &Minerals,ExceptPetroleum
                 'Electronics':'PCU4236--4236--', #HouseholdAppliances& Electrical& ElectronicGoods
                 'Machinery':'PCU4238--4238--', # Machinery,Equipment,& Supplies
                 'Drug':'PCU4242--4242--', # Drugs &Druggists'Sundries
                 'Apparel':'PCU4243--4243--', # Apparel,PieceGoods, & Notions
                 'Grocery':'PCU4244--4244--', # Grocery &RelatedProducts
                 'Chemical':'PCU4246--4246--', # Chemicals& AlliedProducts
                 'Petroleum':'PCU4247--4247--', # Petroleum& PetroleumProducts
                 'Misc':'PCU4249--4249--'
                 }
    
    
    # Year that scalings should be relative to
    year = 2017
    
    output_dir= '/wrk/charkins/emissions/longterm_emis_scaling/Wholesale_trade/20241121/'
    PPI_fn ='/wrk/charkins/emissions/longterm_emis_scaling/PPI/' +'SeriesReport-20241122011736_0d9e4f.xlsx'
    
    pull_production = True
    save_production= True
    pull_PPI = True
    save_PPI = True
    save_adj_production = True

# Calculate production statistics
def calc_production():
    
    # Loop over sales
    sales_data = []
    for commodity in inputs.series_id:
        address= inputs.API_base_address.replace('CATCODE',inputs.series_id[commodity]).replace('DATATYPECODE',inputs.dc_sales)
        print(address)
        r = requests.get(address)
        json_data = r.json()
        df = pd.DataFrame(json_data,#.get('series')[0].get('data'),
                              columns = ['time_slot_id', commodity,'data_type','category_code','error_data','seasonally_adj','time','us'])
        
        df.set_index('time', drop=True, inplace=True)
        df.loc[df.loc[:,commodity] == 'NA',commodity] = np.nan 
        df_out = df.iloc[1:,:]
        df_out = df_out.loc[:,commodity]
        sales_data.append(df_out)
    

    final_sales_data = pd.concat(sales_data, axis=1).fillna(np.nan)
    final_sales_data.sort_index(inplace=True)
    for column in final_sales_data.columns:
        final_sales_data[column]=final_sales_data[column].astype(np.float64)
    print(final_sales_data)
    
    inventory_data = []
    for commodity in inputs.series_id:
        address= inputs.API_base_address.replace('CATCODE',inputs.series_id[commodity]).replace('DATATYPECODE',inputs.dc_inventory)
        print(address)
        r = requests.get(address)
        json_data = r.json()
        df = pd.DataFrame(json_data,#.get('series')[0].get('data'),
                              columns = ['time_slot_id', commodity,'data_type','category_code','error_data','seasonally_adj','time','us'])
        df.set_index('time', drop=True, inplace=True)
        df.loc[df.loc[:,commodity] == 'NA',commodity] = np.nan 
        df_out = df.iloc[1:,:]
        df_out = df_out.loc[:,commodity]
        inventory_data.append(df_out)
    

    final_inventory_data = pd.concat(inventory_data, axis=1).fillna(np.nan)
    final_inventory_data.sort_index(inplace=True)
    for column in final_inventory_data.columns:
        final_inventory_data[column]=final_inventory_data[column].astype(np.float64)
    print(final_inventory_data)
    
    inventory_diff = final_inventory_data.diff()
    
    production = inventory_diff+final_sales_data
    
    print(production)
    
    # Rolling 3 Month Average
    production_rolling = production.rolling(window=3,center=True,axis=0).mean()
    
    # final_data.insert(0,'Year',final_data.index.astype(str).str[:4])
    # final_data.insert(1,'Month',final_data.index.astype(str).str[4:])
    
    production_rolling.insert(0,'Year',production_rolling.index.astype(str).str[:4])
    production_rolling.insert(1,'Month',production_rolling.index.astype(str).str[5:])
    
    print(production_rolling)
    
    print("need to rememeber to adjust sales pr production by the number of days in a month")
    # fuel_scaling = final_data_rolling.copy(deep=True)
    # for fuel in series_id:
    #     fuel_scaling.loc[:,fuel] = final_data_rolling.loc[:,fuel]/(final_data_rolling.loc[final_data_rolling.loc[:,'Year']==str(year),fuel].mean(axis=0))
        
    # for fuel in ['comm','res','egu','ind']:
    #     fuel_scaling.loc[:,fuel] = final_data_rolling.loc[:,fuel]/(final_data_rolling.loc[final_data_rolling.loc[:,'Year']==str(year),fuel].mean(axis=0))
        
    # fuel_scaling.to_csv('fuel_scalings.csv',index=False)
    # final_data_rolling.to_csv('energy_use_3mo_rolling.csv',index=False)
    # final_data.to_csv('energy_use.csv',index=False)
    production.insert(0,'Year',production.index.astype(str).str[:4])
    production.insert(1,'Month',production.index.astype(str).str[5:])
    
    return production
    
# end calc_production

# Read PPI, match up with commodity name and 
def calc_PPI():
    
    # read in the file
    PPI = pd.read_excel(inputs.PPI_fn, sheet_name='BLS Data Series', header=3,
                        dtype={'Series ID':str}).melt(id_vars=['Series ID'],var_name='Date',value_name='producer_price')
    PPI.insert(2,'Year',PPI['Date'].astype(str).str.partition('\n')[2])
    PPI.insert(3,'Month_nam',PPI['Date'].astype(str).str.partition('\n')[0])
    
    
    month_dict = {month: index for index, month in enumerate(calendar.month_abbr) if month}
    months = PPI['Month_nam'].apply(lambda x: month_dict[x])
    PPI.insert(1,'Month',months.astype(str).str.pad(2,fillchar='0'))
    
    # Add in commodity data
    pcu_to_commodity=pd.DataFrame.from_dict(inputs.ppi_id,orient='index',columns=['Series ID']).reset_index().rename(columns={'index':'Commodity'})
    PPI=PPI.join(pcu_to_commodity.set_index('Series ID'),on='Series ID')
    
    # convert all measures to be relative to a given annual average
    for commodity in inputs.ppi_id:
        year_idx= np.logical_and(PPI['Year']==str(inputs.year),PPI['Commodity']==commodity)
        commodity_idx= PPI['Commodity']==commodity
        reference_avg = PPI.loc[year_idx,'producer_price'].mean()
        PPI.loc[commodity_idx,'producer_price'] = PPI.loc[commodity_idx,'producer_price']/reference_avg

    return PPI
# end calc_PPI

# Join PPI and production datasets and calculate adjusted daily production values and 3 month rolling version
def join_datasets(ppi,prod):
    
    # Reformat production and join to PPI
    prod_reformat = prod.melt(id_vars=['Year','Month'],var_name='Commodity',value_name='production')
    df_joined = prod_reformat.join(ppi.set_index(['Commodity','Year','Month'])['producer_price'],on=['Commodity','Year','Month'])
    
    # add column with number of days in each month
    df_joined['num_days'] = df_joined[['Year','Month']].apply(return_num_days,axis=1)
    
    # Calculate adjusted daily production value
    df_joined['prod_adj'] = df_joined['production']/(df_joined['producer_price']*df_joined['num_days'])
    
    # Sort values so we can calculate rolling average
    df_joined.sort_values(by=['Commodity','Year','Month'],inplace=True)
    
    df_joined['prod_adj_3mo_roll']=np.nan
    
    # Loop over commodities to apply rolling
    for commodity in inputs.series_id:
        df_joined.loc[:,'prod_adj_3mo_roll']=df_joined.loc[:,'prod_adj'].rolling(window=3,center=True,axis=0).mean()
        
    
    return df_joined
# end join_datasets

# Returns the number of days in the month
def return_num_days(row):
    return calendar.monthrange(int(row['Year']),int(row['Month']))[1]
# end return_num_days

def main():

    if inputs.pull_production:
        prod=calc_production()
        
    if inputs.pull_PPI:
        ppi = calc_PPI()
        if inputs.save_PPI:
            ppi.to_csv(inputs.output_dir+'PPI.csv',index=False)
    
    df = join_datasets(ppi,prod)
    if inputs.save_adj_production:
        df.to_csv(inputs.output_dir+'wholesale_prod_adjusted.csv',index=False)
    
        

# End Main

if __name__ == "__main__":
    main()


