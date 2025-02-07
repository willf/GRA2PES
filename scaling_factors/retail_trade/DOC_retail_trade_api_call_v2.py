import pandas as pd
import numpy as np
import requests
import calendar

class inputs():
    # https://medium.com/analytics-vidhya/use-python-to-pull-energy-data-from-the-us-department-of-energys-api-11d6f724927e
    
    # Economic indicators API webpage
    # https://www.census.gov/data/developers/data-sets/economic-indicators.html
    # https://www2.census.gov/data/api-documentation/EITS_API_User_Guide_Dec2020.pdf
    
    # Retail trade API address
    # https://api.census.gov/data/timeseries/eits/mrts.html
    
    #XML where you can find variable info links
    
    
    key=''
    # API_base_address = 'https://api.eia.gov/series/?api_key='+key+'&series_id='
    API_base_address = ('https://api.census.gov/data/timeseries/eits/mrts'+
                        '?get=time_slot_id,cell_value&data_type_code=DATATYPECODE&category_code=CATCODE&error_data=no&seasonally_adj=no&for=us:*&time=from+2000&key=' + key)
    
    # Data codes to loop over
    dc_sales = 'SM'
    dc_inventory = 'IM'
    
    series_id = {'Health-Personal':'446', # Health and Personal Care Stores
                 'Pharmacy-DrugStores':'44611',# Pharmacies and Drug Stores 
                 'BuildingMaterials':'4441' # Building Materials and Supplies Dealers
                 } 
    
    # CPI data download from: https://data.bls.gov/cgi-bin/srgate
    # https://www.bls.gov/help/hlpforma.htm
    # https://www.bls.gov/help/hlpforma.htm#CU
    # https://download.bls.gov/pub/time.series/cu/cu.item
    cpi_id= {'Health-Personal':'CUUR0000SA0L1E',
                'Pharmacy-DrugStores':'CUUR0000SA0L1E',} # US City average All items less food and energy }
    # use all items less food+energy
    
    ppi_id = {'BuildingMaterials':'PCU236211236211'}
    
    # Year that scalings should be relative to
    year = 2017
    
    output_dir= '/wrk/charkins/emissions/longterm_emis_scaling/retail_trade/20241121/'
    CPI_fn ='/wrk/charkins/emissions/longterm_emis_scaling/CPI/' +'SeriesReport-20241122012641_fd2366.xlsx'
    PPI_fn ='/wrk/charkins/emissions/longterm_emis_scaling/PPI/' +'SeriesReport-20241122011736_0d9e4f_justconst.xlsx'

    pull_sales=True
    save_sales= True
    pull_CPIPPI=True
    save_CPIPPI = True
    save_adj_sales = True

# Calculate production statistics
def calc_sales():
    
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
    
    final_sales_data.insert(0,'Year',final_sales_data.index.astype(str).str[:4])
    final_sales_data.insert(1,'Month',final_sales_data.index.astype(str).str[5:])
    
    return final_sales_data
    
# end calc_sales

# Read CPI, match up with commodity name and 
def calc_CPI():
    
    # read in the file
    CPI = pd.read_excel(inputs.CPI_fn, sheet_name='BLS Data Series', header=3,
                        dtype={'Series ID':str}).melt(id_vars=['Series ID'],var_name='Date',value_name='consumer_price')
    CPI.insert(2,'Year',CPI['Date'].astype(str).str.partition('\n')[2])
    CPI.insert(3,'Month_nam',CPI['Date'].astype(str).str.partition('\n')[0])
    
    CPI=CPI.loc[np.logical_and(CPI['Month_nam']!='HALF1',CPI['Month_nam']!='HALF2'),:].reset_index(drop=True)
    
    month_dict = {month: index for index, month in enumerate(calendar.month_abbr) if month}
    months = CPI['Month_nam'].apply(lambda x: month_dict[x])
    CPI.insert(1,'Month',months.astype(str).str.pad(2,fillchar='0'))
    
    # Add in commodity data
    pcu_to_commodity=pd.DataFrame.from_dict(inputs.cpi_id,orient='index',columns=['Series ID']).reset_index().rename(columns={'index':'Commodity'})
    CPI=CPI.join(pcu_to_commodity.set_index('Series ID'),on='Series ID')
    
    # convert all measures to be relative to a given annual average
    for commodity in inputs.cpi_id:
        year_idx= np.logical_and(CPI['Year']==str(inputs.year),CPI['Commodity']==commodity)
        commodity_idx= CPI['Commodity']==commodity
        reference_avg = CPI.loc[year_idx,'consumer_price'].mean()
        CPI.loc[commodity_idx,'consumer_price'] = CPI.loc[commodity_idx,'consumer_price']/reference_avg

    return CPI
# end calc_CPI

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
def join_datasets(cpi,prod):
    
    # Reformat production and join to PPI
    prod_reformat = prod.melt(id_vars=['Year','Month'],var_name='Commodity',value_name='sales')
    df_joined = prod_reformat.join(cpi.set_index(['Commodity','Year','Month'])['price'],on=['Commodity','Year','Month'])
    
    # add column with number of days in each month
    df_joined['num_days'] = df_joined[['Year','Month']].apply(return_num_days,axis=1)
    
    # Calculate adjusted daily production value
    df_joined['sales_adj'] = df_joined['sales']/(df_joined['price']*df_joined['num_days'])
    
    # Sort values so we can calculate rolling average
    df_joined.sort_values(by=['Commodity','Year','Month'],inplace=True)
    
    df_joined['sales_adj_3mo_roll']=np.nan
    
    # Loop over commodities to apply rolling
    for commodity in inputs.series_id:
        df_joined.loc[:,'sales_adj_3mo_roll']=df_joined.loc[:,'sales_adj'].rolling(window=3,center=True,axis=0).mean()
        
    
    return df_joined
# end join_datasets

# Returns the number of days in the month
def return_num_days(row):
    return calendar.monthrange(int(row['Year']),int(row['Month']))[1]
# end return_num_days


def main():

    if inputs.pull_sales:
        sales=calc_sales()
        
    if inputs.pull_CPIPPI:
        cpi = calc_CPI()
        ppi = calc_PPI()
        cpi= cpi.rename(columns={'consumer_price':'price'})
        ppi= ppi.rename(columns={'producer_price':'price'})
        cpippi = pd.concat([cpi,ppi])
        if inputs.save_CPIPPI:
            cpippi.to_csv(inputs.output_dir+'CPI-PPI.csv',index=False)
    print(sales)
    
    df = join_datasets(cpippi,sales)
    if inputs.save_adj_sales:
        df.to_csv(inputs.output_dir+'retail_sales_adjusted.csv',index=False)
    
        

# End Main

if __name__ == "__main__":
    main()


