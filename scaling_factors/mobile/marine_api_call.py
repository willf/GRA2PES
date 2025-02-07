import pandas as pd
import numpy as np
import requests
import calendar

class inputs():
    # https://www.census.gov/foreign-trade/reference/guides/Guide%20to%20International%20Trade%20Datasets.pdf
    # https://www.census.gov/foreign-trade/reference/definitions/index.html#SWT
    
    # https://api.census.gov/data/timeseries/intltrade/exports/statehs/variables.html
    
    # Variables
    # STATE - 2-character State of Origin of Movement, using US Postal Service State Abbreviations. "XX" = Unidentified
    # YEAR - 4-character Year
    # MONTH - 	2-character Month
    # VES_WGT_MO_FLAG - 15-digit Vessel Shipping Weight, missing flag
    # VES_WGT_MO - 15-digit Vessel Shipping Weight
    
    colnames = ['STATE','YEAR','MONTH','VES_WGT_MO','time']
    
    # types 
    trade_type=['exports','imports']
    
    API_base_address = ('https://api.census.gov/data/timeseries/intltrade/[TRADETYPE]/statehs?get=STATE,YEAR,MONTH,VES_WGT_MO&time=from+2005-12+to+2024-03')
    
    output_dir= '/wrk/charkins/emissions/longterm_emis_scaling/mobile/20241121/'

    pull_ship=True
    save_ship=True
    
#end inputs

# Calculate production statistics
def calc_ship():
    
    data = []
    
    for trade in inputs.trade_type:
        address= inputs.API_base_address.replace('[TRADETYPE]',trade)
        print(address)
        r = requests.get(address)
        json_data = r.json()
        df = pd.DataFrame(json_data,#.get('series')[0].get('data'),
                              columns = inputs.colnames)
        df.rename(columns={'YEAR':'Year','MONTH':'Month','VES_WGT_MO':trade},inplace=True)
        df = df[1:]
        df[trade]=df[trade].astype(np.float64)
        data.append(df)
        
    all_data = data[0]
    all_data=all_data.join(data[1].set_index(['time','STATE'])['imports'],on=['time','STATE'])
    all_data=all_data.loc[all_data['STATE']!='-',:]
    
    final_data = all_data.groupby(['Year','Month','time'])[['imports','exports']].sum().reset_index()
    final_data=final_data.sort_values(by=['Year','Month'])
    final_data['num_days'] = final_data[['Year','Month']].apply(return_num_days,axis=1)
    final_data['ship_total_adj']= (final_data['imports']+final_data['exports'])/final_data['num_days']
    final_data['ship_total_3mo_roll']=final_data['ship_total_adj'].rolling(window=3,center=True,axis=0).mean()
    
    return final_data
    
# end calc_rail

# Returns the number of days in the month
def return_num_days(row):
    return calendar.monthrange(int(row['Year']),int(row['Month']))[1]
# end return_num_days

def main():

    if inputs.pull_ship:
        ship_traffic=calc_ship()
        
    if inputs.save_ship:
        ship_traffic.to_csv(inputs.output_dir+'ship_traffic_adjusted.csv',index=False)
    
        

# End Main

if __name__ == "__main__":
    main()
