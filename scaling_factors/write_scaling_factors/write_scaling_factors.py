# Code to generate annual pollutant scaling factors from NEI tier summaries 

import pandas as pd
import numpy as np

class inputs():
    #input_dir = 'C:/Users/CH/OneDrive/research/longterm_emis/write_scaling_factors/'
    #output_dir_area = 'C:/Users/CH/OneDrive/research/longterm_emis/write_scaling_factors/AREA_scaling_factors/'
    #output_dir_pt = 'C:/Users/CH/OneDrive/research/longterm_emis/write_scaling_factors/POINT_scaling_factors/'
    
    input_dir = '/wrk/charkins/emissions/longterm_emis_scaling/write_scaling_factors/'
    output_dir_area = '/wrk/charkins/emissions/longterm_emis_scaling/write_scaling_factors/AREA_scaling_factors_202411/'
    output_dir_pt = '/wrk/charkins/emissions/longterm_emis_scaling/write_scaling_factors/POINT_scaling_factors_202411/'
    
    # information for mobile scalings (rail,ship, aviation)
    mobile_scaling_dir = '/wrk/charkins/emissions/longterm_emis_scaling/mobile/20241121/'
    rail_data_fn = mobile_scaling_dir + 'rail_traffic_adjusted.csv'
    plane_data_fn = mobile_scaling_dir + 'aviation_traffic_adjusted.csv'
    ship_data_fn =  mobile_scaling_dir + 'ship_traffic_adjusted.csv'
    
    # Information for wholesale trade scalings
    wholesale_scaling_dir = '/wrk/charkins/emissions/longterm_emis_scaling/Wholesale_trade/20241121/'
    wholesale_data_fn = wholesale_scaling_dir + 'wholesale_prod_adjusted.csv'
    
    # Information for retail trade scalings
    retail_scaling_dir = '/wrk/charkins/emissions/longterm_emis_scaling/retail_trade/20241121/'
    retail_data_fn = retail_scaling_dir + 'retail_sales_adjusted.csv'
    
    # Information for fuel consumption scalings
    energy_scaling_dir = '/wrk/charkins/emissions/longterm_emis_scaling/Energy/20241121/'
    energy_data_fn = energy_scaling_dir + 'energy_use_3mo_rolling.csv'

    
    # This is the year the scaling should be relative to
    origin_year = 2017
    origin_year_VCP_McD = 2012
     
    # sector files to make 
    area_sectors_out = ['COMM','FUG','INDF','INDP','MOB','OG','RES',
                        'VCP','VCP_McD']
    point_sectors_out = ['COMM','EGU','FUG','INDF','INDP','MOB','OG',
                        'VCP']
    
    write_years = [2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023]
    
    time_vars = ['Year','Month']
    
    write_area = True
    save_area = True
    area_fn_out = 'Area[SECTOR].csv'
    write_point = True
    save_point = True
    point_fn_out = 'Pt[SECTOR].csv'
    
# end inputs

# Combines all of the dataframes into one dataframe
def combine_all_df(inputs,rail_df,ship_df,plane_df,wholesale_df,
               retail_df,energy_df):
    
    # Construct the output dataframe to join data to
    df_out = pd.DataFrame({'Year':sorted(inputs.write_years*12),
                           'Month':[1,2,3,4,5,6,7,8,9,10,11,12]*len(inputs.write_years)})
    
    time_vars = inputs.time_vars

    # Join rail data
    rail_vars = ['rail_total_3mo_roll'] 
    df_out = df_out.join(rail_df[time_vars + rail_vars].set_index(time_vars),on=time_vars)
    
    # Join ship data
    ship_vars = ['ship_total_3mo_roll'] 
    df_out = df_out.join(ship_df[time_vars + ship_vars].set_index(time_vars),on=time_vars)
    
    # Join plane data
    plane_vars = ['aviation_3mo_roll'] 
    df_out = df_out.join(plane_df[time_vars + plane_vars].set_index(time_vars),on=time_vars)
    
    # Join wholesale data
    wholesale_vars = ['TotalMerchant','MotorVehicle','Lumber','ProfessionalSupplies',
                      'Metals','Electronics','Machinery','Drug','Apparel',
                      'Grocery','Chemical','Petroleum','Misc'] 
    # First need to reshape the wholesale data so we have columns for each data type
    wholesale_df = wholesale_df[['Year','Month','Commodity','prod_adj_3mo_roll']].set_index(['Year','Month']).pivot(columns='Commodity')['prod_adj_3mo_roll'].reset_index()
    df_out = df_out.join(wholesale_df[time_vars + wholesale_vars].set_index(time_vars),on=time_vars)
    
    # Join retail data
    retail_vars = ['BuildingMaterials','Health-Personal','Pharmacy-DrugStores'] 
    # First need to reshape the wholesale data so we have columns for each data type
    retail_df = retail_df[['Year','Month','Commodity','sales_adj_3mo_roll']].set_index(['Year','Month']).pivot(columns='Commodity')['sales_adj_3mo_roll'].reset_index()
    df_out = df_out.join(retail_df[time_vars + retail_vars].set_index(time_vars),on=time_vars)
    
    # Join energy data
    energy_vars = ['res','comm','res_NG','comm_NG',
                   'ind_wood','ind_coal','ind_NG','ind_distillate',
                   'egu_wood','egu_coal','egu_NG','egu_distillate'] 
    df_out = df_out.join(energy_df[time_vars + energy_vars].set_index(time_vars),on=time_vars)
    
    return df_out
# end

# Constructs all of the datasets for each scaling group, then saves to .csv
def make_area(inputs,df):
    
    for sector in inputs.area_sectors_out:
        make_area_subsectors(inputs,sector=sector,df=df)

# end make_area

# Contains code to make scalings for each sector 
def make_area_subsectors(inputs,sector,df):

    time_vars = inputs.time_vars
    base = df[time_vars].copy(deep=True)    

    if sector=='COMM':
        df_out = base
        comm = df['comm']/(df.loc[df['Year']==inputs.origin_year,'comm'].mean())
        df_out['COMM_Coal'] = comm
        df_out['COMM_Oil'] = comm
        df_out['COMM_NG'] = comm
        df_out['COMM_Wood'] = comm
    elif sector=='FUG': 
        df_out = base
        petro = df['Petroleum']/(df.loc[df['Year']==inputs.origin_year,'Petroleum'].mean())
        lumber = df['Lumber']/(df.loc[df['Year']==inputs.origin_year,'Lumber'].mean())
        df_out['OffG'] = petro
        df_out['STORAGE'] = petro
        df_out['EVAPGAS'] = petro
        df_out['CONST'] = lumber
        df_out['DUST'] = lumber
    elif sector=='INDF': 
        df_out = base
        df_out['IND_Coal'] = df['ind_coal']/(df.loc[df['Year']==inputs.origin_year,'ind_coal'].mean())
        df_out['IND_Oil'] = df['ind_distillate']/(df.loc[df['Year']==inputs.origin_year,'ind_distillate'].mean())
        df_out['IND_NG'] = df['ind_NG']/(df.loc[df['Year']==inputs.origin_year,'ind_NG'].mean())
        df_out['IND_Wood'] = df['ind_wood']/(df.loc[df['Year']==inputs.origin_year,'ind_wood'].mean())
    elif sector=='INDP': 
        df_out = base
        df_out['CHEM'] = df['Chemical']/(df.loc[df['Year']==inputs.origin_year,'Chemical'].mean())
        df_out['FOOD'] = df['Grocery']/(df.loc[df['Year']==inputs.origin_year,'Grocery'].mean())
        df_out['METAL'] = df['Metals']/(df.loc[df['Year']==inputs.origin_year,'Metals'].mean())
        df_out['PULP'] = df['Lumber']/(df.loc[df['Year']==inputs.origin_year,'Lumber'].mean())
        df_out['MACHINE'] = df['Machinery']/(df.loc[df['Year']==inputs.origin_year,'Machinery'].mean())
        df_out['OTH'] = df['Misc']/(df.loc[df['Year']==inputs.origin_year,'Misc'].mean())
        df_out['MISC'] = df['Misc']/(df.loc[df['Year']==inputs.origin_year,'Misc'].mean())
    elif sector=='MOB': 
        df_out = base
        ship = df['ship_total_3mo_roll']/(df.loc[df['Year']==inputs.origin_year,'ship_total_3mo_roll'].mean())
        df_out['CMV_DSL'] = ship
        df_out['CMV_OGV'] = ship
        df_out['RAIL'] = df['rail_total_3mo_roll']/(df.loc[df['Year']==inputs.origin_year,'rail_total_3mo_roll'].mean())
    elif sector=='OG': 
        df_out = base
        petro = df['Petroleum']/(df.loc[df['Year']==inputs.origin_year,'Petroleum'].mean())
        NG = (df['ind_NG']+df['res_NG']+ df['egu_NG']+df['comm_NG'])/((
            df.loc[df['Year']==inputs.origin_year,'ind_NG']+
            df.loc[df['Year']==inputs.origin_year,'res_NG']+
            df.loc[df['Year']==inputs.origin_year,'egu_NG']+
            df.loc[df['Year']==inputs.origin_year,'comm_NG']).mean())
        OnG = (petro+NG)/2 # Take the average scaling of petroleum and NG consumption
        df_out['OnG'] = OnG
    elif sector=='RES':
        df_out = base
        res = df['res']/(df.loc[df['Year']==inputs.origin_year,'res'].mean())
        df_out['RES_Coal'] = res
        df_out['RES_Oil'] = res
        df_out['RES_NG'] = res
        df_out['RES_Wood'] = res
    elif sector=='VCP':
        df_out = base
        chem = df['Chemical']/(df.loc[df['Year']==inputs.origin_year,'Chemical'].mean())
        df_out['IndCoat'] = chem
        df_out['Degreasing'] = chem
        df_out['Inks'] = chem
        df_out['IndAdhesive'] = chem
        df_out['AgPesticide'] = np.ones(chem.shape)
        df_out['TotlVCP'] = chem
    elif sector=='VCP_McD':
        df_out = base
        chem = df['Chemical']/(df.loc[df['Year']==inputs.origin_year_VCP_McD,'Chemical'].mean())
        const = df['BuildingMaterials']/(df.loc[df['Year']==inputs.origin_year_VCP_McD,'BuildingMaterials'].mean())
        personal = (df['Health-Personal']-df['Pharmacy-DrugStores'])/((
            df.loc[df['Year']==inputs.origin_year_VCP_McD,'Health-Personal']-
            df.loc[df['Year']==inputs.origin_year_VCP_McD,'Pharmacy-DrugStores']).mean())
        df_out['ConsCoatings'] = const
        df_out['IndCoatings'] = chem
        df_out['Inks'] = chem
        df_out['ConsAdhesives'] = const
        df_out['IndAdhesives'] = chem
        df_out['Personal'] = personal
        df_out['Cleaning'] = personal

    # save file
    if inputs.save_area:
        fn_out = inputs.output_dir_area + inputs.area_fn_out.replace('[SECTOR]',sector)
        df_out.to_csv(fn_out,index=False)
    
# end make_area_subsectors

# Constructs all of the datasets for each scaling group, then saves to .csv
def make_point(inputs,df):
    
    for sector in inputs.point_sectors_out:
        make_point_subsectors(inputs,sector=sector,df=df)

# end make_area

# Contains code to make scalings for each sector 
def make_point_subsectors(inputs,sector,df):

    time_vars = inputs.time_vars
    base = df[time_vars].copy(deep=True)    

    if sector=='COMM':
        df_out = base
        comm = df['comm']/(df.loc[df['Year']==inputs.origin_year,'comm'].mean())
        df_out['PtCOMM_Coal'] = comm
        df_out['PtCOMM_Oil'] = comm
        df_out['PtCOMM_NG'] = comm
        df_out['PtCOMM_Wood'] = comm
    elif sector=='EGU':
        df_out = base
        df_out['PtEGU_Coal'] = df['egu_coal']/(df.loc[df['Year']==inputs.origin_year,'egu_coal'].mean())
        df_out['PtEGU_NG'] = df['egu_NG']/(df.loc[df['Year']==inputs.origin_year,'egu_NG'].mean())
        df_out['PtEGU_Oil'] = df['egu_distillate']/(df.loc[df['Year']==inputs.origin_year,'egu_distillate'].mean())
        df_out['PtEGU_BIO'] = df['egu_wood']/(df.loc[df['Year']==inputs.origin_year,'egu_wood'].mean())
    elif sector=='FUG': 
        df_out = base
        petro = df['Petroleum']/(df.loc[df['Year']==inputs.origin_year,'Petroleum'].mean())
        lumber = df['Lumber']/(df.loc[df['Year']==inputs.origin_year,'Lumber'].mean())
        df_out['PtSTORAGE'] = petro
        df_out['PtEVAPGAS'] = petro
        df_out['PtCONST'] = lumber
    elif sector=='INDF': 
        df_out = base
        df_out['PtIND_Coal'] = df['ind_coal']/(df.loc[df['Year']==inputs.origin_year,'ind_coal'].mean())
        df_out['PtIND_NG'] = df['ind_NG']/(df.loc[df['Year']==inputs.origin_year,'ind_NG'].mean())
        df_out['PtIND_NG2'] = df['ind_NG']/(df.loc[df['Year']==inputs.origin_year,'ind_NG'].mean())
        df_out['PtIND_Oil'] = df['ind_distillate']/(df.loc[df['Year']==inputs.origin_year,'ind_distillate'].mean())
        df_out['PtIND_Oil2'] = df['ind_distillate']/(df.loc[df['Year']==inputs.origin_year,'ind_distillate'].mean())
        df_out['PtIND_BIO'] = df['ind_wood']/(df.loc[df['Year']==inputs.origin_year,'ind_wood'].mean())
    elif sector=='INDP': 
        df_out = base
        df_out['PtCHEM'] = df['Chemical']/(df.loc[df['Year']==inputs.origin_year,'Chemical'].mean())
        df_out['PtFOOD'] = df['Grocery']/(df.loc[df['Year']==inputs.origin_year,'Grocery'].mean())
        df_out['PtMETAL'] = df['Metals']/(df.loc[df['Year']==inputs.origin_year,'Metals'].mean())
        df_out['PtREFINE'] = df['Petroleum']/(df.loc[df['Year']==inputs.origin_year,'Petroleum'].mean())
        df_out['PtPULP'] = df['Lumber']/(df.loc[df['Year']==inputs.origin_year,'Lumber'].mean())
        df_out['PtELECT'] = df['Electronics']/(df.loc[df['Year']==inputs.origin_year,'Electronics'].mean())
        df_out['PtMOTOR'] = df['MotorVehicle']/(df.loc[df['Year']==inputs.origin_year,'MotorVehicle'].mean())
        df_out['PtAPPAREL'] = df['Apparel']/(df.loc[df['Year']==inputs.origin_year,'Apparel'].mean())
        df_out['PtPHOTO'] = df['ProfessionalSupplies']/(df.loc[df['Year']==inputs.origin_year,'ProfessionalSupplies'].mean())
        df_out['PtDRUG'] = df['Drug']/(df.loc[df['Year']==inputs.origin_year,'Drug'].mean())
        df_out['PtMISC'] = df['TotalMerchant']/(df.loc[df['Year']==inputs.origin_year,'TotalMerchant'].mean())
        df_out['PtMISC2'] = df['TotalMerchant']/(df.loc[df['Year']==inputs.origin_year,'TotalMerchant'].mean())
    elif sector=='MOB': 
        df_out = base
        df_out['PtAVIATION'] = df['aviation_3mo_roll']/(df.loc[df['Year']==inputs.origin_year,'aviation_3mo_roll'].mean())
        df_out['PtRAIL'] = df['rail_total_3mo_roll']/(df.loc[df['Year']==inputs.origin_year,'rail_total_3mo_roll'].mean())
    elif sector=='OG': 
        df_out = base
        petro = df['Petroleum']/(df.loc[df['Year']==inputs.origin_year,'Petroleum'].mean())
        NG = (df['ind_NG']+df['res_NG']+ df['egu_NG']+df['comm_NG'])/((
            df.loc[df['Year']==inputs.origin_year,'ind_NG']+
            df.loc[df['Year']==inputs.origin_year,'res_NG']+
            df.loc[df['Year']==inputs.origin_year,'egu_NG']+
            df.loc[df['Year']==inputs.origin_year,'comm_NG']).mean())
        OnG = (petro+NG)/2 # Take the average scaling of petroleum and NG consumption
        df_out['PtOnG'] = OnG
    elif sector=='VCP':
        df_out = base
        chem = df['Chemical']/(df.loc[df['Year']==inputs.origin_year,'Chemical'].mean())
        df_out['PtAdhesives'] = chem
        df_out['PtCoatings'] = chem
        df_out['PtDegreasing'] = chem
        df_out['PtInks'] = chem
        
    # save file
    if inputs.save_point:
        fn_out = inputs.output_dir_pt + inputs.point_fn_out.replace('[SECTOR]',sector)
        df_out.to_csv(fn_out,index=False)
    
# end make_point_subsectors

def main():
    # Open datasets that will be used
    rail_df = pd.read_csv(inputs.rail_data_fn)
    ship_df = pd.read_csv(inputs.ship_data_fn)
    plane_df = pd.read_csv(inputs.plane_data_fn)
    wholesale_df = pd.read_csv(inputs.wholesale_data_fn)
    retail_df = pd.read_csv(inputs.retail_data_fn)
    energy_df = pd.read_csv(inputs.energy_data_fn)
    
    # Combine all of the columns so it's easy to access everything
    combined_df = combine_all_df(inputs,rail_df,ship_df,plane_df,wholesale_df,
                   retail_df,energy_df)
    
    if inputs.write_area:
        make_area(inputs,combined_df)
    if inputs.write_point:
        make_point(inputs,combined_df)
    
# End Main

if __name__ == "__main__":
    main()
