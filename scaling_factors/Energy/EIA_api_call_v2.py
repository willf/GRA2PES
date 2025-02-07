import pandas as pd
import numpy as np
import requests

def main():
    # https://medium.com/analytics-vidhya/use-python-to-pull-energy-data-from-the-us-department-of-energys-api-11d6f724927e
    # API Key
    key = ''
    
    # Energy Consumption by Sector
    # https://www.eia.gov/opendata/qb.php?category=711226
    
    #API_base_address = 'https://api.eia.gov/series/?api_key='+key+'&series_id='
    
    #API_base_address = 'https://api.eia.gov/v2/total-energy/data/?frequency=monthly&data[0]=value&facets[msn][]=BMRCBUS&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000'
    API_base_address = 'https://api.eia.gov/v2/total-energy/data/?frequency=monthly&data[0]=value&facets[msn][]=$SERIES_ID&sort[0][column]=period&sort[0][direction]=desc&api_key='+key
    
    #series_id = {'res_wood':'TOTAL.BMRCBUS.M','res_NG':'TOTAL.NNRCBUS.M','res_distillate':'TOTAL.PARCBUS.M',# Residential Monthly Energy (Trillion Btu)
    #             'comm_wood':'TOTAL.BMCCBUS.M','comm_coal':'TOTAL.CLCCBUS.M','comm_NG':'TOTAL.NNCCBUS.M','comm_distillate':'TOTAL.PMCCBUS.M', # Commercial Monthly Energy (Trillion Btu)
    #             'ind_wood':'TOTAL.BMICBUS.M','ind_coal':'TOTAL.CLICBUS.M','ind_NG':'TOTAL.NNICBUS.M','ind_distillate':'TOTAL.PMICBUS.M', # Industrial Monthly Energy (Trillion Btu)
    #             'egu_wood':'TOTAL.BMEIBUS.M','egu_coal':'TOTAL.CLEIBUS.M','egu_NG':'TOTAL.NNEIBUS.M','egu_distillate':'TOTAL.PAEIBUS.M'} # Power Generation Monthly Energy (Trillion Btu)
                 
    series_id = {'res_wood':'BMRCBUS','res_NG':'NNRCBUS','res_distillate':'PARCBUS',# Residential Monthly Energy (Trillion Btu)
                 'comm_wood':'BMCCBUS','comm_coal':'CLCCBUS','comm_NG':'NNCCBUS','comm_distillate':'PMCCBUS', # Commercial Monthly Energy (Trillion Btu)
                 'ind_wood':'BMICBUS','ind_coal':'CLICBUS','ind_NG':'NNICBUS','ind_distillate':'PMICBUS', # Industrial Monthly Energy (Trillion Btu)
                 'egu_wood':'BMEIBUS','egu_coal':'CLEIBUS','egu_NG':'NNEIBUS','egu_distillate':'PAEIBUS'} # Power Generation Monthly Energy (Trillion Btu)
    # # Year that scalings should be relative to
    # year = 2017
    
    output_dir = '/wrk/charkins/emissions/longterm_emis_scaling/Energy/20241121/'
    
    
    data = []
    for fuel in series_id:
        r = requests.get(API_base_address.replace("$SERIES_ID",series_id[fuel]))
        print(API_base_address.replace("$SERIES_ID",series_id[fuel]))
        json_data = r.json()
        #df = pd.DataFrame(json_data.get('series')[0].get('data'),
        #                      columns = ['Date', fuel])
        df = pd.DataFrame(json_data.get('response')['data'],
                              columns = ['period', 'value']).rename(columns=
                              {'period':'Date','value':fuel}).astype({fuel:'float'})
        
        df.set_index('Date', drop=True, inplace=True)
        df.loc[df.loc[:,fuel] == 'NA',fuel] = np.nan 
        data.append(df)
    
    final_data = pd.concat(data, axis=1).fillna(np.nan)
    final_data.sort_index(inplace=True)
    
    final_data['res'] = final_data['res_distillate'] + final_data['res_wood'] + + final_data['res_NG'] 
    final_data['comm'] = final_data['comm_distillate'] + final_data['comm_coal'] + final_data['comm_wood'] + final_data['comm_NG'] 
    final_data['egu'] = final_data['egu_distillate'] + final_data['egu_coal'] + final_data['egu_wood'] + final_data['egu_NG'] 
    final_data['ind'] = final_data['ind_distillate'] + final_data['ind_coal'] + final_data['ind_wood'] + final_data['ind_NG'] 
    
    # Rolling 3 Month Average
    final_data_rolling = final_data.rolling(window=3,center=True,axis=0).mean()
    
    final_data.insert(0,'Year',final_data.index.astype(str).str[:4])
    final_data.insert(1,'Month',final_data.index.astype(str).str[4:])
    
    final_data_rolling.insert(0,'Year',final_data_rolling.index.astype(str).str[:4])
    final_data_rolling.insert(1,'Month',final_data_rolling.index.astype(str).str[5:])
    #breakpoint()
    #final_data.insert(0,'Year',final_data.index.astype(str).str.split('-')[0])
    #final_data.insert(1,'Month',final_data.index.astype(str).str.split('-')[1])
    
    #final_data_rolling.insert(0,'Year',final_data_rolling.index.astype(str).str.split('-')[0])
    #final_data_rolling.insert(1,'Month',final_data_rolling.index.astype(str).str.split('-')[1])
    
    # fuel_scaling = final_data_rolling.copy(deep=True)
    # for fuel in series_id:
    #     fuel_scaling.loc[:,fuel] = final_data_rolling.loc[:,fuel]/(final_data_rolling.loc[final_data_rolling.loc[:,'Year']==str(year),fuel].mean(axis=0))
        
    # for fuel in ['comm','res','egu','ind']:
    #     fuel_scaling.loc[:,fuel] = final_data_rolling.loc[:,fuel]/(final_data_rolling.loc[final_data_rolling.loc[:,'Year']==str(year),fuel].mean(axis=0))
        
    # fuel_scaling.to_csv('fuel_scalings.csv',index=False)
    final_data.to_csv(output_dir+'energy_use.csv',index=False)
    final_data_rolling.to_csv(output_dir+'energy_use_3mo_rolling.csv',index=False)
    # final_data.to_csv('energy_use.csv',index=False)

# End Main

if __name__ == "__main__":
    main()


