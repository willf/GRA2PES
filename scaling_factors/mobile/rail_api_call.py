import pandas as pd
import numpy as np
import requests
import calendar

class inputs():
    # https://data.bts.gov/Research-and-Statistics/Monthly-Transportation-Statistics/crem-w557
    # https://dev.socrata.com/foundry/data.bts.gov/crem-w557
    
    
    
    
    
    
    # key='d5bc9e5f9eff425eb5ddb810ce6de87a07b1a50f'
    # API_base_address = 'https://api.eia.gov/series/?api_key='+key+'&series_id='
    API_base_address = ('https://data.bts.gov/resource/crem-w557.json?$select=date, system_use_freight_rail_1, system_use_freight_rail')
    
    output_dir= '/wrk/charkins/emissions/longterm_emis_scaling/mobile/20241121/'

    pull_rail=True
    save_rail=True

# Calculate production statistics
def calc_rail():
    

    address= inputs.API_base_address
    print(address)
    r = requests.get(address)
    json_data = r.json()
    df = pd.DataFrame(json_data,#.get('series')[0].get('data'),
                          columns = ['date', 'system_use_freight_rail_1','system_use_freight_rail'])
    df.rename(columns={'system_use_freight_rail_1':'intermodal_units','system_use_freight_rail':'carloads'},inplace=True)
    
    df['carloads'] = df['carloads'].astype(np.float64)
    df['intermodal_units'] = df['intermodal_units'].astype(np.float64)
    
    df.insert(0,'Year',df['date'].astype(str).str[:4])
    df.insert(1,'Month',df['date'].astype(str).str[5:7])
    df['num_days'] = df[['Year','Month']].apply(return_num_days,axis=1)
    df['rail_total']= (df['intermodal_units']+df['carloads'])/df['num_days']
    
    df['rail_total_3mo_roll']=df['rail_total'].rolling(window=3,center=True,axis=0).mean()
    
    return df
    
# end calc_rail




# Returns the number of days in the month
def return_num_days(row):
    return calendar.monthrange(int(row['Year']),int(row['Month']))[1]
# end return_num_days


def main():

    if inputs.pull_rail:
        rail_traffic=calc_rail()
        
    if inputs.save_rail:
        rail_traffic.to_csv(inputs.output_dir+'rail_traffic_adjusted.csv',index=False)
    
        

# End Main

if __name__ == "__main__":
    main()


