# Uses geopandas_env
import numpy as np
import pandas as pd
# import geopandas as gpd
# import shapely
import calendar
import time

class inputs():
    input_dir = '/wrk/d2/charkins/FOG_data/Production/'
    output_dir = '/wrk/d2/charkins/FOG_data/Production/Production_formatted/'
    
    
    state_list = ['AL','AR','AZ','CA','CO','FL',
                  'KS','KY','LA','MD','MI','MO',
                  'MS','MT','ND','NE','NM','NV',
                  'NY','OH','OK','OR','PA','SD',
                  'TN','TX','UT','VA','WV','WY']
    # input filenames in input_dir
    raw_data_fn = '[ST] Producing Entity Monthly Production.CSV'
    well_data_fn = '[ST] Production Table.CSV'
    
    file_out = '[ST]_Monthly_Production_YYYY_Month##.csv'
    subtotals_file_out = '[ST]_Production_totals_YYYY.csv'
    
    # Create file for what time period
    monthly=True # if monthly is true look for specific months data
    months=[1,2,3,4,5,6,7,8,9,10,11,12]
    # months=[1,2,3,4,5,6,7,8,9,10,11,12]
    years=[2022,2023]
    
    DyCO2_ratio = 0.65748591500000
    LcCO2_ratio=  1.17499154700000

    save_output = True
    
# end inputs class

def main():
    
    start_time = time.time()
    print('START OF PROGRAM: ')
    print('--------------------------------------------------')
    
    state_totals = pd.DataFrame(columns=['State','Year','Month','Oprod','Gprod']) 
    
    for state in inputs.state_list:
        start_state_time = time.time()
        print('Formatting state production files for State: ' +state)
        
        # Read in production and well info and join all
        prod=pd.read_csv(inputs.input_dir+inputs.raw_data_fn.replace('[ST]',state),
                         usecols=['Entity ID','Monthly Production Date','Monthly Oil', 'Monthly Gas'],
                         dtype={'Entity ID':str,'Monthly Production Date':str,'Monthly Oil':np.float64, 'Monthly Gas':np.float64})
        prod=prod.rename(columns={'Monthly Oil':'Oprod','Monthly Gas':'Gprod'})
        
        well_info=pd.read_csv(inputs.input_dir+inputs.well_data_fn.replace('[ST]',state),
                              usecols=['Entity ID','Surface Latitude (WGS84)', 'Surface Longitude (WGS84)'],
                              dtype={'Entity ID':str,'Surface Latitude (WGS84)':np.float64, 'Surface Longitude (WGS84)':np.float64})
        well_info=well_info.rename(columns={'Surface Latitude (WGS84)':'Lat','Surface Longitude (WGS84)':'Lon'})
        
        prod['Year']=pd.DatetimeIndex(prod['Monthly Production Date']).year
        prod['Month']=pd.DatetimeIndex(prod['Monthly Production Date']).month
        for year in inputs.years:
            
            # if monthly return only production from the given month
            if inputs.monthly:
                totals = state_totals.copy(deep=True)
                for month in inputs.months:
                    print('Working on files for: ' + str(year) + ', Month : ' + str(month) + ', State: ' + state )
                    
                    num_days = calendar.monthrange(year,month)[1]
                    prod_out = prod.loc[np.logical_and(prod['Month']==month,prod['Year']==year),:].fillna(0).reset_index(drop=True) 
                    # Convert to daily production
                    prod_out['Oprod']=prod_out['Oprod']/num_days
                    prod_out['Gprod']=prod_out['Gprod']/num_days
                    prod_matched=prod_out.join(well_info.set_index('Entity ID'),on='Entity ID')
                    prod_matched['DyCO2'] = prod_matched['Gprod']*inputs.DyCO2_ratio
                    prod_matched['LcCO2'] = prod_matched['Gprod']*inputs.LcCO2_ratio

                    if inputs.save_output:
                        fn = inputs.output_dir+inputs.file_out.replace('[ST]',state).replace('YYYY',str(year)).replace('##',str(month).zfill(2))
                        prod_matched[['Lat','Lon','Oprod','Gprod','DyCO2','LcCO2']].to_csv(fn,index=False)
                    
                    Oprod_totl = prod_out['Oprod'].sum()
                    Gprod_totl = prod_out['Gprod'].sum()
                    totals = totals.append({'State':state, 'Year':year,
                                            'Month':month, 'Oprod':Oprod_totl, 
                                            'Gprod':Gprod_totl }, ignore_index=True)
                #end for over months
                # save out monthly production totals
                totals = totals.append({'State':state, 'Year':year,
                                        'Month':0, 'Oprod':totals['Oprod'].sum()/12, 
                                        'Gprod':totals['Gprod'].sum()/12 },ignore_index=True)
                fn_totals = inputs.output_dir+inputs.subtotals_file_out.replace('[ST]',state).replace('YYYY',str(year))
                totals.to_csv(totals.to_csv(fn_totals,index=False))
                
                
                    
                    
            # if not monthly, return total production
            else:
                num_days=365
                prod_out=prod.fillna(0).groupby(['Entity ID']).sum().reset_index()
                prod_matched=prod_out.join(well_info.set_index('Entity ID'),on='Entity ID')
                prod_out['Oprod']=prod_out['Oprod']/num_days
                prod_out['Gprod']=prod_out['Gprod']/num_days
                prod_matched['DyCO2'] = prod_matched['Gprod']*inputs.DyCO2_ratio
                prod_matched['LcCO2'] = prod_matched['Gprod']*inputs.LcCO2_ratio
                
                if inputs.save_output:
                    fn = inputs.output_dir+inputs.file_out.replace('[ST]',state).replace('YYYY',str(year)).replace('##','00')
                    prod_matched[['Lat','Lon','Oprod','Gprod','DyCO2','LcCO2']].to_csv(fn,index=False)
            # end if monthly
        # end loop over years
        end_state_time = time.time()
        # print('Time to make files for Year : ' + str(year) + ',  Month: ' + str(m) + ' is ' + str( round(end_state_time-start_state_time ) )  + ' s')
        print('Time elapsed to make files for Year : ' + str(year) + ', State: ' + state + ' is ' + str( round(end_state_time-start_state_time ) )  + ' s')
    # end loop over states
    print('END OF PROGRAM: ')
    print('--------------------------------------------------')
    end_time = time.time()
    print('Total elapsed time (s): ' + str(round(end_time-start_time)))        
    

if __name__ == "__main__":
    main()
