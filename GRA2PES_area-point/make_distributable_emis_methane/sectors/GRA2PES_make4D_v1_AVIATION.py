# Uses geopandas_env
import xarray as xr
import numpy as np
import pandas as pd
# import datetime as dt
from netCDF4 import Dataset
from wrf import getvar, interpline, CoordPair, xy_to_ll, ll_to_xy
import time
import gc

class inputs():
    # Input/Output Dirs
    
    gridded_files = []#'/wrk/users/charkins/emissions/V7_GRA2PES/AREA[YY]_ncf/Month[MM]/AreaRAIL/[DOW]/AreaRAIL_[dayhalf].nc']
    
    # CANNOT HAVE MORE THAN ONE CURRENTLY 
    point_files = ['/wrk/users/charkins/emissions/V7_GRA2PES/POINT[YY]_ncf/Month[MM]/PtAVIATION/[DOW]/PtAVIATION_[dayhalf].nc']
    
    out_file1 = '/wrk/users/charkins/emissions/V7_GRA2PES/GRA2PESv1.0_methane/GRA2PESv1.0_AVIATION/[YYYY][MM]/[DOW]/GRA2PESv1.0_AVIATION_[YYYY][MM]_[DOW]_00to11Z.nc'
    out_file2 = '/wrk/users/charkins/emissions/V7_GRA2PES/GRA2PESv1.0_methane/GRA2PESv1.0_AVIATION/[YYYY][MM]/[DOW]/GRA2PESv1.0_AVIATION_[YYYY][MM]_[DOW]_12to23Z.nc'
    
    wrfchemi_example = '/wrk/charkins/emissions/total_emis/4km_domain/wrfinput_d01_conus4k'
    
    skip_methane_output = False
    
    sector_title_output = 'Aviation Emissions'
    
    years = [2021]
    months = [1,2,3,4,5,6,7,8,9,10,11,12]
    dows = ['weekdy','satdy','sundy']
    
class out_nc_file():
    
    # unit conversions to get into final units
    vars_unit_conv = {    'ffCO2':10**6/44.01/16,
                'CO2':10**6/44.01/16,
                'CO':10**6/28.01/16,
                'NH3':10**6/17.031/16,
                'NOX':10**6/46.0055/16,
                'SO2':10**6/64.066/16,
                'VOC':1/16,
                'PM25-PRI':1/16,
                'PM10-PRI':1/16,
                'HC01':1/16,
                'HC02':1/16,
                'HC03':1/16,
                'HC04':1/16,
                'HC05':1/16,
                'HC06':1/16,
                'HC07':1/16,
                'HC08':1/16,
                'HC09':1/16,
                'HC10':1/16,
                'HC11':1/16,
                'HC12':1/16,
                'HC13':1/16,
                'HC14':1/16,
                'HC15':1/16,
                'HC16':1/16,
                'HC17':1/16,
                'HC18':1/16,
                'HC19':1/16,
                'HC20':1/16,
                'HC21':1/16,
                'HC22':1/16,
                'HC23':1/16,
                'HC24':1/16,
                'HC25':1/16,
                'HC26':1/16,
                'HC27':1/16,
                'HC28':1/16,
                'HC29':1/16,
                'HC30':1/16,
                'HC31':1/16,
                'HC32':1/16,
                'HC33':1/16,
                'HC34':1/16,
                'HC35':1/16,
                'HC36':1/16,
                'HC37':1/16,
                'HC38':1/16,
                'HC39':1/16,
                'HC40':1/16,
                'HC41':1/16,
                'HC42':1/16,
                'HC43':1/16,
                'HC44':1/16,
                'HC45':1/16,
                'HC46':1/16,
                'HC47':1/16,
                'HC48':1/16,
                'HC49':1/16,
                'HC50':1/16,
                'HC51':1/16,
                'HC52':1/16,
                'HC53':1/16,
                'HC54':1/16,
                'HC55':1/16,
                'HC56':1/16,
                'HC57':1/16,
                'HC58':1/16,
                'HC59':1/16,
                'HC60':1/16,
                'HC61':1/16,
                'HC62':1/16,
                'HC63':1/16,
                'HC64':1/16,
                'HC65':1/16,
                'HC66':1/16,
                'HC67':1/16,
                'HC68':1/16,
                'PM01':1/16,
                'PM02':1/16,
                'PM03':1/16,
                'PM04':1/16,
                'PM05':1/16,
                'PM06':1/16,
                'PM07':1/16,
                'PM08':1/16,
                'PM09':1/16,
                'PM10':1/16,
                'PM11':1/16,
                'PM12':1/16,
                'PM13':1/16,
                'PM14':1/16,
                'PM15':1/16,
                'PM16':1/16,
                'PM17':1/16,
                'PM18':1/16,
                'PM19':1/16,
                }
    
    # variable description info
    variables = {    
                #'ffCO2':{'units':'mole km^-2 hr^-1','description':'Fossil Fuel Carbon Dioxide'},
                # 'CO2':{'units':'mole km^-2 hr^-1','description':'Carbon Dioxide'},
                # 'CO':{'units':'mole km^-2 hr^-1','description':'Carbon Monoxide'},
                # 'NH3':{'units':'mole km^-2 hr^-1','description':'Ammonia'},
                # 'NOX':{'units':'mole km^-2 hr^-1','description':'Nitrogen Oxides (NOX)'},
                # 'SO2':{'units':'mole km^-2 hr^-1','description':'Sulfur Dioxide'},
                # 'VOC':{'units':'metric tons km^-2 hr^-1','description':'Total NMVOC Emissions'},
                # 'PM25-PRI':{'units':'metric tons km^-2 hr^-1','description':'Total Primary PM2.5 Emissions'},
                # 'PM10-PRI':{'units':'metric tons km^-2 hr^-1','description':'Total Primary PM10 Emissions'},
                'HC01':{'units':'mole km^-2 hr^-1','description':'Methane'},
                # 'HC02':{'units':'mole km^-2 hr^-1','description':'Ethane+Alkanes with k(OH)< 500 ppm^-1-min^-1'},
                # 'HC03':{'units':'mole km^-2 hr^-1','description':'Alkanes with 500<k(OH)<2500 ppm^-1-min^-1,EXCLUDE(C3H8,C2H2,ethanol,acids)'},
                # 'HC04':{'units':'mole km^-2 hr^-1','description':'Alkanes with 2500<k(OH)<5000 ppm^-1-min^-1,EXCLUDE(butanes)'},
                # 'HC05':{'units':'mole km^-2 hr^-1','description':'Alkanes with 5000<k(OH)<10000 ppm^-1-min^-1,EXCLUDE(pentanes)'},
                # 'HC06':{'units':'mole km^-2 hr^-1','description':'Alkanes with k(OH)>10000 ppm^-1-min^-1,EXCLUDE(ethylene glycol)'},
                # 'HC07':{'units':'mole km^-2 hr^-1','description':'Ethylene'},
                # 'HC08':{'units':'mole km^-2 hr^-1','description':'Alkenes with k(OH)<20000 ppm^-1-min^-1'},
                # 'HC09':{'units':'mole km^-2 hr^-1','description':'Alkenes with k(OH)>20000 ppm^-1-min^-1,EXCLUDE(dienes and styrenes)'},
                # 'HC10':{'units':'mole km^-2 hr^-1','description':'Anthropogenic Isoprene'},
                # 'HC11':{'units':'mole km^-2 hr^-1','description':'Anthropogenic Terpenes (VCPs)'},
                # 'HC12':{'units':'mole km^-2 hr^-1','description':'Aromatics with k(OH)<20000 ppm^-1-min^-1,EXCLUDE(benzene and toluene)'},
                # 'HC13':{'units':'mole km^-2 hr^-1','description':'Aromatics with k(OH)<20000 ppm^-1-min^-1,EXCLUDE((xylenes)'},
                # 'HC14':{'units':'mole km^-2 hr^-1','description':'Formaldehyde'},
                # 'HC15':{'units':'mole km^-2 hr^-1','description':'Acetaldehyde'},
                # 'HC16':{'units':'mole km^-2 hr^-1','description':'>C2 aldehydes'},
                # 'HC17':{'units':'mole km^-2 hr^-1','description':'Benzaldehyde'},
                # 'HC18':{'units':'mole km^-2 hr^-1','description':'Acetone'},
                # 'HC19':{'units':'mole km^-2 hr^-1','description':'Methylethyl keton'},
                # 'HC20':{'units':'mole km^-2 hr^-1','description':'PRD2 SAPRAC species (aromatic ketones)'},
                # 'HC21':{'units':'mole km^-2 hr^-1','description':'Methanol'},
                # 'HC22':{'units':'mole km^-2 hr^-1','description':'Glyoxal'},
                # 'HC23':{'units':'mole km^-2 hr^-1','description':'Methylglyoxal'},
                # 'HC24':{'units':'mole km^-2 hr^-1','description':'Biacetyl'},
                # 'HC25':{'units':'mole km^-2 hr^-1','description':'Phenols'},
                # 'HC26':{'units':'mole km^-2 hr^-1','description':'Cresols'},
                # 'HC27':{'units':'mole km^-2 hr^-1','description':'Methacrolein'},
                # 'HC28':{'units':'mole km^-2 hr^-1','description':'Methylvinyl ketone'},
                # 'HC29':{'units':'mole km^-2 hr^-1','description':'IPRD SAPRAC species (>C4 unsaturated aldehydes)'},
                # 'HC30':{'units':'mole km^-2 hr^-1','description':'Formic Acid'},
                # 'HC31':{'units':'mole km^-2 hr^-1','description':'Acetic Acid'},
                # 'HC32':{'units':'mole km^-2 hr^-1','description':'>C2 Acids'},
                # 'HC33':{'units':'mole km^-2 hr^-1','description':'Xylenols'},
                # 'HC34':{'units':'mole km^-2 hr^-1','description':'Catechols'},
                # 'HC35':{'units':'mole km^-2 hr^-1','description':'NonVolatile Compounds'},
                # 'HC36':{'units':'mole km^-2 hr^-1','description':'Propylene'},
                # 'HC37':{'units':'mole km^-2 hr^-1','description':'Acetylene'},
                # 'HC38':{'units':'mole km^-2 hr^-1','description':'Benzene'},
                # 'HC39':{'units':'mole km^-2 hr^-1','description':'Butanes'},
                # 'HC40':{'units':'mole km^-2 hr^-1','description':'Pentanes'},
                # 'HC41':{'units':'mole km^-2 hr^-1','description':'Toluene'},
                # 'HC42':{'units':'mole km^-2 hr^-1','description':'m-Xylene'},
                # 'HC43':{'units':'mole km^-2 hr^-1','description':'o-Xylene'},
                # 'HC44':{'units':'mole km^-2 hr^-1','description':'p-Xylene'},
                # 'HC45':{'units':'mole km^-2 hr^-1','description':'Propane'},
                # 'HC46':{'units':'mole km^-2 hr^-1','description':'Dienes'},
                # 'HC47':{'units':'mole km^-2 hr^-1','description':'Styrenes'},
                # 'HC48':{'units':'mole km^-2 hr^-1','description':'Ethanol'},
                # 'HC49':{'units':'mole km^-2 hr^-1','description':'Ethylene Glycol'},
                # 'HC50':{'units':'mole km^-2 hr^-1','description':'Unidentified_Unknown VOC'},
                # 'HC51':{'units':'mole km^-2 hr^-1','description':'Isopropyl Alcohol (OVCP)'},
                # 'HC52':{'units':'mole km^-2 hr^-1','description':'Propylene Glycol (OVCP)'},
                # 'HC53':{'units':'mole km^-2 hr^-1','description':'Glycerol (OVCP)'},
                # 'HC54':{'units':'mole km^-2 hr^-1','description':'D4-Siloxane'},
                # 'HC55':{'units':'mole km^-2 hr^-1','description':'D5-Siloxane'},
                # 'HC56':{'units':'mole km^-2 hr^-1','description':'Other Siloxane'},
                # 'HC57':{'units':'mole km^-2 hr^-1','description':'NROG'},
                # 'HC58':{'units':'mole km^-2 hr^-1','description':'PCBTF'},
                # 'HC59':{'units':'mole km^-2 hr^-1','description':'PDCBZ'},
                # 'HC60':{'units':'mole km^-2 hr^-1','description':'Propanal'},
                # 'HC61':{'units':'mole km^-2 hr^-1','description':'Butanal'},
                # 'HC62':{'units':'mole km^-2 hr^-1','description':'Pentanal'},
                # 'HC63':{'units':'mole km^-2 hr^-1','description':'Hexanal'},
                # 'HC64':{'units':'mole km^-2 hr^-1','description':'Heptanal'},
                # 'HC65':{'units':'mole km^-2 hr^-1','description':'Octanal'},
                # 'HC66':{'units':'mole km^-2 hr^-1','description':'Nonanal'},
                # 'HC67':{'units':'mole km^-2 hr^-1','description':'Unsaturated Aldehydes'},
                # 'HC68':{'units':'mole km^-2 hr^-1','description':'C10+ Aldehydes'},
                # 'PM01':{'units':'metric tons km^-2 hr^-1','description':'Nonspeciated Primary PM2.5 (sum of PM08-PM19)'},
                # 'PM02':{'units':'metric tons km^-2 hr^-1','description':'Sulfate PM2.5'},
                # 'PM03':{'units':'metric tons km^-2 hr^-1','description':'Nitrate PM2.5'},
                # 'PM04':{'units':'metric tons km^-2 hr^-1','description':'Organic Carbon PM2.5'},
                # 'PM05':{'units':'metric tons km^-2 hr^-1','description':'Elemental Carbon PM2.5'},
                # 'PM06':{'units':'metric tons km^-2 hr^-1','description':'Non-Carbon Organic PM2.5'},
                # 'PM07':{'units':'metric tons km^-2 hr^-1','description':'Ammonium PM2.5'},
                # 'PM08':{'units':'metric tons km^-2 hr^-1','description':'Aluminum PM2.5'},
                # 'PM09':{'units':'metric tons km^-2 hr^-1','description':'Calcium PM2.5'},
                # 'PM10':{'units':'metric tons km^-2 hr^-1','description':'Iron PM2.5'},
                # 'PM11':{'units':'metric tons km^-2 hr^-1','description':'Water PM2.5'},
                # 'PM12':{'units':'metric tons km^-2 hr^-1','description':'Magnesium PM2.5'},
                # 'PM13':{'units':'metric tons km^-2 hr^-1','description':'Other PM2.5'},
                # 'PM14':{'units':'metric tons km^-2 hr^-1','description':'Potassium PM2.5'},
                # 'PM15':{'units':'metric tons km^-2 hr^-1','description':'Manganese PM2.5'},
                # 'PM16':{'units':'metric tons km^-2 hr^-1','description':'Chloride PM2.5'},
                # 'PM17':{'units':'metric tons km^-2 hr^-1','description':'Sodium PM2.5'},
                # 'PM18':{'units':'metric tons km^-2 hr^-1','description':'Titanium PM2.5'},
                # 'PM19':{'units':'metric tons km^-2 hr^-1','description':'Silicon PM2.5'},
                }

# end inputs class

class data():
    vert_grid_edges = [0,16.8,50.5,84.3,127,170,256,343,476,610,791,883,975,1160,1350,1644,1943,2252,2677,3010,3350] # interval grid level height (m) of wind-data (climatology) used in momentum lift calcs
    
    #wind speed at level height (m) specified by REFWZ, first dim is time, second dim is level
    wind_speed = np.array([
      [ 4.27, 4.01, 4.16, 4.27, 4.30, 4.26, 4.20, 4.12, 4.11, 4.11, 4.15, 4.25, 4.44, 3.43, 3.61, 3.86, 4.12, 4.36, 4.62, 4.94], 
      [ 5.12, 5.21, 5.24, 5.12, 5.30, 5.62, 5.79, 5.96, 6.02, 5.96, 5.87, 5.74, 5.73, 5.70, 5.75, 5.85, 5.96, 4.51, 4.67, 4.93], 
      [ 5.26, 5.58, 5.95, 6.40, 6.70, 6.88, 6.99, 6.91, 6.14, 6.82, 7.07, 7.27, 7.36, 7.28, 7.19, 7.02, 6.99, 6.92, 6.93, 7.02], 
      [ 7.08, 5.46, 5.43, 5.64, 5.97, 6.31, 6.72, 7.24, 7.60, 7.88, 8.09, 8.12, 6.87, 7.73, 8.08, 8.36, 8.45, 8.37, 8.28, 8.09], 
      [ 8.03, 7.92, 7.87, 7.91, 7.91, 6.34, 6.15, 6.22, 6.50, 6.82, 7.25, 7.80, 8.22, 8.58, 8.88, 9.04, 7.53, 8.46, 8.87, 9.22], 
      [ 9.34, 9.26, 9.16, 8.95, 8.86, 8.70, 8.60, 8.57, 8.51, 6.94, 6.75, 6.74, 6.94, 7.21, 7.61, 8.16, 8.62, 9.05, 9.43, 9.73], 
      [ 8.17, 9.20, 9.72,10.14,10.30,10.26,10.16, 9.92, 9.77, 9.55, 9.36, 9.25, 9.11, 7.42, 7.28, 7.38, 7.52, 7.65, 7.95, 8.44], 
      [ 8.92, 9.42, 9.92,10.38, 8.55, 9.65,10.28,10.81,11.03,11.02,10.92,10.65,10.44,10.17, 9.93, 9.77, 9.59, 7.72, 7.72, 7.95], 
      [ 8.29, 8.37, 8.35, 8.64, 9.07, 9.57,10.12,10.71, 8.74, 9.82,10.49,11.08,11.36,11.39,11.32,11.05,10.83,10.54,10.25,10.12], 
      [10.02, 7.89, 7.96, 8.21, 8.69, 8.89, 8.78, 8.88, 9.24, 9.71,10.26,10.89, 8.82, 9.84,10.50,11.11,11.41,11.43,11.33,11.02], 
      [10.75,10.44,10.18,10.14,10.25, 7.92, 8.09, 8.39, 8.93, 9.32, 9.34, 9.32, 9.59, 9.96,10.45,11.06, 8.86, 9.81,10.42,10.98], 
      [11.20,11.11,10.92,10.54,10.23, 9.91, 9.71, 9.77,10.02, 7.95, 8.20, 8.62, 9.23, 9.72, 9.93, 9.92,10.03,10.31,10.70,11.19], 
      [ 8.89, 9.74,10.28,10.76,10.88,10.72,10.46,10.07, 9.80, 9.55, 9.39, 9.50, 9.81, 8.09, 8.36, 8.83, 9.48,10.02,10.33,10.34], 
      [10.40,10.60,10.90,11.26, 8.90, 9.67,10.14,10.54,10.61,10.42,10.16, 9.79, 9.56, 9.34, 9.24, 9.40, 9.73, 8.19, 8.45, 8.94], 
      [ 9.63,10.21,10.56,10.58,10.65,10.79,11.02,11.31, 8.96, 9.56, 9.91,10.21,10.18, 9.96, 9.71, 9.39, 9.25, 9.12, 9.12, 9.33], 
      [ 9.69, 8.38, 8.63, 9.17, 9.91,10.54,10.91,11.01,11.03,11.12,11.22,11.40, 9.10, 9.50, 9.66, 9.84, 9.77, 9.58, 9.42, 9.23], 
      [ 9.27, 9.31, 9.36, 9.53, 9.83, 8.71, 8.98, 9.53,10.30,10.93,11.33,11.47,11.45,11.48,11.51,11.52, 9.45, 9.64, 9.60, 9.62], 
      [ 9.54, 9.44, 9.42, 9.39, 9.58, 9.74, 9.79, 9.89,10.08, 9.24, 9.49, 9.96,10.70,11.36,11.77,11.96,11.95,11.90,11.83,11.68], 
      [10.05,10.08, 9.90, 9.80, 9.72, 9.74, 9.90,10.02,10.28,10.45,10.47,10.52,10.63,10.05,10.23,10.54,11.15,11.76,12.14,12.28], 
      [12.32,12.30,12.21,11.97,10.66,10.66,10.50,10.37,10.32,10.42,10.66,10.81,11.05,11.22,11.25,11.29,11.35,10.99,11.08,11.20], 
      [11.62,12.10,12.37,12.44,12.52,12.62,12.62,12.43,11.40,11.53,11.45,11.42,11.43,11.50,11.68,11.77,11.98,12.21,12.33,12.40], 
      [12.46,12.17,12.18,12.08,12.23,12.44,12.54,12.61,12.84,13.13,13.33,13.34,12.21,12.47,12.43,12.47,12.51,12.53,12.63,12.70], 
      [12.95,13.26,13.44,13.56,13.64,13.37,13.36,13.07,12.91,12.82,12.71,12.74,13.07,13.52,13.93,14.19,12.76,13.18,13.21,13.32], 
      [13.40,13.43,13.54,13.64,13.97,14.32,14.53,14.68,14.72,14.37,14.37,14.05,13.78,13.56,13.31,13.17,13.42,13.87,14.33,14.73]])   

def point_coords_to_xy(lat,lon):
    
    ncfile = Dataset(inputs.wrfchemi_example)
    
    x_y = ll_to_xy(wrfin=ncfile,latitude=lat,longitude=lon)
    
    return x_y

def calc_stack_level(stack_height):
    vert_grid_max = max(data.vert_grid_edges)
    stack_level_func = lambda h,grid_max: next(i for i,v in enumerate(data.vert_grid_edges) if v > h)-1 if h < grid_max else len(data.vert_grid_edges)-2
    return xr.apply_ufunc(stack_level_func,stack_height,vert_grid_max,vectorize=True)

def calc_dhm(stack_diameter, stack_velocity,stack_level,hour):
    dhm_func = lambda d,v,l,h: 3*d*v/data.wind_speed[h,l]
    return xr.apply_ufunc(dhm_func,stack_diameter, stack_velocity,stack_level,hour,vectorize=True)

def calc_fracs(wrfchemi_ds,file):
    
    point_ds = xr.open_dataset(file).fillna({'STKht':0,'STKve':1,'STKdiam':1})
    
    
    x_y = point_coords_to_xy(lat=point_ds['XLAT'].values,lon=point_ds['XLONG'].values).rename({'idx':'ROW'})
    
    
    stack_level=calc_stack_level(point_ds['STKht']) # stack level on vertical grid
    
    dhm1 = xr.zeros_like(point_ds['HC01'],dtype=float)
    dhm1.attrs={}
    dhm2 = xr.zeros_like(point_ds['HC01'],dtype=float)
    dhm2.attrs={}
    
    for hour in np.arange(0,12):
        dhm1.loc[dict(Time=hour)] = calc_dhm(point_ds['STKdiam'],point_ds['STKve'],stack_level,hour)
    for hour in np.arange(12,24):
        dhm2.loc[dict(Time=hour-12)] = calc_dhm(point_ds['STKdiam'],point_ds['STKve'],stack_level,hour)
    
    top1=point_ds['STKht']+1.5*dhm1
    top2=point_ds['STKht']+1.5*dhm2
    
    bot1=point_ds['STKht']+0.5*dhm1
    bot2=point_ds['STKht']+0.5*dhm2
    
    kbot1 = xr.zeros_like(bot1,dtype=float)
    kbot2 = xr.zeros_like(bot2,dtype=float)
    
    ktop1 = xr.zeros_like(top1,dtype=float)
    ktop2 = xr.zeros_like(top2,dtype=float)
    
    
    
    # Create empty array with correct shape to hold vertical info
    vert_rows = xr.concat([xr.zeros_like(point_ds['STKht']).expand_dims(dim="bottom_top_stag")]*len(data.vert_grid_edges),dim="bottom_top_stag")
    vert_inds = xr.concat([xr.ones_like(point_ds['STKht']).expand_dims(dim="bottom_top_stag")]*len(data.vert_grid_edges),dim="bottom_top_stag")
    
    for i in np.arange(len(data.vert_grid_edges)):
        vert_inds.loc[dict(bottom_top_stag=i)]=vert_inds.loc[dict(bottom_top_stag=i)]*i
    
    
    for r in np.arange(0,kbot1.sizes['ROW']):
        x = x_y.loc[dict(x_y='x',ROW=r)].data
        y = x_y.loc[dict(x_y='y',ROW=r)].data
        # Check if in domain
        if 0<=x<wrfchemi_ds.sizes['west_east'] and 0<=y<wrfchemi_ds.sizes['south_north']:
            vert_info = wrfchemi_ds['h_agl'].squeeze().loc[dict(south_north=y,west_east=x,bottom_top_stag=slice(0,len(data.vert_grid_edges)))]
            vert_rows.loc[dict(ROW=r)]=vert_info
            
    
    
    kbot1 = xr.where(vert_rows > bot1, vert_inds-1,len(data.vert_grid_edges)-2).min(dim=['bottom_top_stag']) # set to max possible value so min doesn't catch it unless bot1 is above the max of vert_rows
    ktop1 = xr.where(vert_rows > top1, vert_inds-1,len(data.vert_grid_edges)-2).min(dim=['bottom_top_stag'])
    ktop1 = xr.where(kbot1 >= ktop1, kbot1+1, ktop1)
    ktop1 = xr.where(ktop1 > len(data.vert_grid_edges)-2, len(data.vert_grid_edges)-2, ktop1) # fix possibility of dim out of range created by last line and just put all in highest level if needed
    
    kbot2 = xr.where(vert_rows > bot2, vert_inds-1,len(data.vert_grid_edges)-2).min(dim=['bottom_top_stag']) # set to max possible value so min doesn't catch it unless bot1 is above the max of vert_rows
    ktop2 = xr.where(vert_rows > top2, vert_inds-1,len(data.vert_grid_edges)-2).min(dim=['bottom_top_stag'])
    ktop2 = xr.where(kbot2 >= ktop2, kbot2+1, ktop2)
    ktop2 = xr.where(ktop2 > len(data.vert_grid_edges)-2, len(data.vert_grid_edges)-2, ktop2) # fix possibility of dim out of range created by last line and just put all in highest level if needed
    
    zdif1 = xr.zeros_like(kbot1,dtype=float)
    zdif2 = xr.zeros_like(kbot2,dtype=float)
    
    for hour in np.arange(0,12):
        ind_bottom_top_stag_kbot1 = xr.DataArray(kbot1.loc[dict(Time=hour)].data, dims=['d']).astype(int)
        ind_row_kbot1 = xr.DataArray(np.arange(0,kbot1.sizes['ROW']), dims=['d']).astype(int)
        ind_bottom_top_stag_ktop1 = xr.DataArray(ktop1.loc[dict(Time=hour)].data+1, dims=['d']).astype(int)
        ind_row_ktop1 = xr.DataArray(np.arange(0,ktop1.sizes['ROW']), dims=['d']).astype(int)
        
        ind_bottom_top_stag_kbot2 = xr.DataArray(kbot2.loc[dict(Time=hour)].data, dims=['d']).astype(int)
        ind_row_kbot2 = xr.DataArray(np.arange(0,kbot2.sizes['ROW']), dims=['d']).astype(int)
        ind_bottom_top_stag_ktop2 = xr.DataArray(ktop2.loc[dict(Time=hour)].data+1, dims=['d']).astype(int)
        ind_row_ktop2 = xr.DataArray(np.arange(0,ktop2.sizes['ROW']), dims=['d']).astype(int)
        
        zdif1.loc[dict(Time=hour)] =  vert_rows[ind_bottom_top_stag_ktop1,ind_row_ktop1].rename({'d':'ROW'}) - vert_rows[ind_bottom_top_stag_kbot1,ind_row_kbot1].rename({'d':'ROW'})
        zdif2.loc[dict(Time=hour)] =  vert_rows[ind_bottom_top_stag_ktop2,ind_row_ktop2].rename({'d':'ROW'}) - vert_rows[ind_bottom_top_stag_kbot2,ind_row_kbot2].rename({'d':'ROW'})
    
    
    vert_rows_diff = vert_rows.diff(dim="bottom_top_stag").rename({'bottom_top_stag':'bottom_top'})
    
    # initialize creation of fraction in each grid
    frac_sub1 = xr.ones_like(vert_rows_diff,dtype=float)
    frac_sub2 = xr.ones_like(vert_rows_diff,dtype=float)
    for i in np.arange(len(data.vert_grid_edges)-1):
        frac_sub1.loc[dict(bottom_top=i)]=frac_sub1.loc[dict(bottom_top=i)]*i
        frac_sub2.loc[dict(bottom_top=i)]=frac_sub2.loc[dict(bottom_top=i)]*i
    frac_places1 = xr.concat([frac_sub1.expand_dims("Time")]*12,dim="Time")
    frac_places1 = xr.where(frac_places1>=kbot1,frac_places1,-1)
    frac_places1 = xr.where(frac_places1<=ktop1,frac_places1,-1)
    frac_places1 = xr.where(frac_places1>=0,1,0)
    frac1 = frac_places1*vert_rows_diff/zdif1
    
    frac_places2 = xr.concat([frac_sub2.expand_dims("Time")]*12,dim="Time")
    frac_places2 = xr.where(frac_places2>=kbot2,frac_places2,-1)
    frac_places2 = xr.where(frac_places2<=ktop2,frac_places2,-1)
    frac_places2 = xr.where(frac_places2>=0,1,0)
    frac2 = frac_places2*vert_rows_diff/zdif2

    return frac1,frac2, x_y

def make_output_ds_structure(wrfchemi_ds):
    
    out = wrfchemi_ds[['XLAT','XLONG','PH','PHB','h_agl','h_asl']].squeeze().loc[dict(bottom_top_stag=slice(0,21))].copy(deep=True)
    
    var_structure = xr.concat([xr.zeros_like(wrfchemi_ds['QKE'],dtype=float)]*12,dim='Time').loc[dict(bottom_top=slice(0,20))]
    
    out['example_variable']=xr.zeros_like(var_structure,dtype=float)
    encoding_ph = {'dtype': 'float32', 'chunksizes':(wrfchemi_ds.sizes['south_north'], wrfchemi_ds.sizes['west_east']),
                      'zlib': True, 'complevel': 1, '_FillValue': None,'coordinates':'XLONG XLAT'}
    encoding_phb = {'dtype': 'float32', 'chunksizes':(1,wrfchemi_ds.sizes['south_north'], wrfchemi_ds.sizes['west_east']),
                      'zlib': True, 'complevel': 1, '_FillValue': None,'coordinates':'XLONG XLAT'}
    out['PH'].encoding=encoding_ph
    out['PHB'].encoding=encoding_phb
    
    for v in ['XLAT','XLONG','PH','PHB']:
        del out[v].attrs['FieldType']
        del out[v].attrs['MemoryOrder']
    
    out['h_asl'].attrs={'description': 'Height above sea level', 'units': 'm', 'stagger': 'Z'}
    out['h_agl'].attrs={'description': 'Height above ground level', 'units': 'm', 'stagger': 'Z'}
    
    # NEED TO SET GLOBAL ATTRIBUTES
    extract_keys = {'TITLE','DX','DY','GRID_ID','PARENT_ID','CEN_LAT','CEN_LON','TRUELAT1','TRUELAT2','MOAD_CEN_LAT','STAND_LON','POLE_LAT',
        'POLE_LON','MAP_PROJ','MAP_PROJ_CHAR'}
    out.attrs={key: out.attrs[key] for key in out.attrs.keys()& extract_keys}
    out.attrs['TITLE']='GRA2PESv1.0 Emissions Dataset'
    
    return out
    
def sum_files_variable(files,year,month,dow,dayhalf,variable):
    
    if variable == 'ffCO2':
        drop_var_list = [*(set([*out_nc_file.variables.keys()]) - set([variable,'CO2']))] + ['Times','time']
    else:
        drop_var_list = [*(set([*out_nc_file.variables.keys()]) - set([variable]))] + ['Times','time']
    
    for c,file in enumerate(files):
        fn = file.replace('[YYYY]',str(year)).replace('[YY]',str(year)[2:]).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow).replace('[dayhalf]',dayhalf)
        
        if c == 0: 
            ds_out = xr.open_dataset(fn,cache=False,drop_variables=drop_var_list)
            if 'ffCO2' not in ds_out.data_vars and variable == 'ffCO2':
                ds_out['ffCO2'] = ds_out['CO2'].copy(deep=True)
            
            if 'XLAT' in ds_out.data_vars:
                ds_out = ds_out.set_coords(['XLAT'])
            if 'XLONG' in ds_out.data_vars:
                ds_out = ds_out.set_coords(['XLONG'])
            
            
            if not ds_out.data_vars:
                if 'ROW' in ds_out.dims:
                    ds_out[variable] =  xr.DataArray(data=np.zeros((12,ds_out.sizes['ROW'])),dims=['ROW'])
                else:
                    ds_out[variable] = xr.DataArray(data=np.zeros((12,ds_out.sizes['south_north'],ds_out.sizes['west_east'])),dims=['Time','south_north','west_east'])
            
            
           
        else:
            ds = xr.open_dataset(fn,cache=False,drop_variables=drop_var_list)
            if 'ffCO2' not in ds.data_vars and variable == 'ffCO2':
                ds['ffCO2'] = ds['CO2'].copy(deep=True)
            
            if 'XLAT' in ds.data_vars:
                ds = ds.set_coords(['XLAT'])
            if 'XLONG' in ds.data_vars:
                ds = ds.set_coords(['XLONG'])
            
            if not ds.data_vars:
                if 'ROW' in ds.dims:
                    ds[variable] =  xr.DataArray(data=np.zeros((12,ds.sizes['ROW'])),dims=['ROW'])
                else:
                    ds[variable] = xr.DataArray(data=np.zeros((12,ds.sizes['south_north'],ds.sizes['west_east'])),dims=['Time','south_north','west_east'])
                
            
            
            ds_out[variable] = ds_out[variable] + ds[variable]
            ds.close()
    
    return ds_out
    
def place_vars(year,month,dow,out_ds_empty,frac1,frac2,x_y):
    
    
    for file in inputs.gridded_files + inputs.point_files:
        fn1 = file.replace('[YYYY]',str(year)).replace('[YY]',str(year)[2:]).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow).replace('[dayhalf]','00to12Z')
        fn2 = file.replace('[YYYY]',str(year)).replace('[YY]',str(year)[2:]).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow).replace('[dayhalf]','12to24Z')
        print('Summing file: ' + fn1)
        print('Summing file: ' + fn2)
    
    
    for c, variable in enumerate(out_nc_file.variables.keys()):
        
        # Check if methane output should be skipped for now
        if variable == 'HC01' and inputs.skip_methane_output:
            continue
        
        print('Working on emis for variable: '+ variable )
        time_loop_variable_start = time.time()
        
        if inputs.gridded_files:
            print('Making gridded dataset 00z: ')
            area_emis1 = sum_files_variable(inputs.gridded_files,year=year,month=month,dow=dow,dayhalf='00to12Z',variable=variable)
            print('Making gridded dataset 12z: ')
            area_emis2 = sum_files_variable(inputs.gridded_files,year=year,month=month,dow=dow,dayhalf='12to24Z',variable=variable)
        
        if inputs.point_files:
            print('Making point dataset 00z: ')
            point_emis1 = sum_files_variable(inputs.point_files,year=year,month=month,dow=dow,dayhalf='00to12Z',variable=variable).assign_coords({'x':x_y.loc[dict(x_y='x')].drop_vars(['x_y','latlon_coord']),'y':x_y.loc[dict(x_y='y')].drop_vars(['x_y','latlon_coord'])})
            print('Making point dataset 12z: ')
            point_emis2 = sum_files_variable(inputs.point_files,year=year,month=month,dow=dow,dayhalf='12to24Z',variable=variable).assign_coords({'x':x_y.loc[dict(x_y='x')].drop_vars(['x_y','latlon_coord']),'y':x_y.loc[dict(x_y='y')].drop_vars(['x_y','latlon_coord'])})
        
        
        if inputs.point_files:
            point_emis_sum_variables = point_emis1.sum(dim=['Time','ROW']) + point_emis2.sum(dim=['Time','ROW'])
        if inputs.gridded_files:
            area_emis_sum_variables = area_emis1.sum(dim=['Time','south_north','west_east']) + area_emis2.sum(dim=['Time','south_north','west_east'])
        
        
        emis1= out_ds_empty.copy(deep=True).rename({'example_variable':variable})
        emis2 = out_ds_empty.copy(deep=True).rename({'example_variable':variable})
        
        
        # Make the time variable 
        date = str(year)+'-'+str(month).zfill(2)+'-'+'01'
        emis1.coords['Time'] = pd.date_range(date+" 00:00:00",date+" 12:00:00",freq='h',inclusive='left')
        emis2.coords['Time'] = pd.date_range(date+" 12:00:00",date+" 23:01:00",freq='h',inclusive='left')
        
        hrs_00z = ['_00:00:00','_01:00:00','_02:00:00','_03:00:00','_04:00:00','_05:00:00',
                   '_06:00:00','_07:00:00','_08:00:00','_09:00:00','_10:00:00','_11:00:00']
        hrs_12z = ['_12:00:00','_13:00:00','_14:00:00','_15:00:00','_16:00:00','_17:00:00',
                   '_18:00:00','_19:00:00','_20:00:00','_21:00:00','_22:00:00','_23:00:00']
        dates_00z = [i+j for i,j in zip([date]*12,hrs_00z)]
        dates_12z = [i+j for i,j in zip([date]*12,hrs_12z)]
        emis1['Times'] = xr.DataArray(data=dates_00z,dims=['Time']).astype('unicode')
        emis1['Times'].encoding={'char_dim_name':'DateStrLen'}
        emis2['Times'] = xr.DataArray(data=dates_12z,dims=['Time']).astype('unicode')
        emis2['Times'].encoding={'char_dim_name':'DateStrLen'}
        
        unit_conv = out_nc_file.vars_unit_conv[variable]
        
        if inputs.point_files:
            if variable in point_emis1.data_vars:
                
                # Do the stuff for point sources
                if point_emis_sum_variables[variable]>0:
                    emis_points1_grouped = aggregate_by_xy(point_emis1[variable]*frac1*unit_conv,variable,out_ds_empty.sizes)
                    emis_points2_grouped = aggregate_by_xy(point_emis2[variable]*frac2*unit_conv,variable,out_ds_empty.sizes)
                    
                    ind_x1 = xr.DataArray((emis_points1_grouped['x'].values.astype(int)).tolist(),dims=["d"])
                    ind_y1 = xr.DataArray((emis_points1_grouped['y'].values.astype(int)).tolist(),dims=["d"])
                    ind_t1 = xr.DataArray((emis_points1_grouped['Time'].values.astype(int)).tolist(),dims=["d"])
                    ind_bottom_top1 = xr.DataArray((emis_points1_grouped['bottom_top'].values.astype(int)).tolist(),dims=["d"])
                    
                    ind_x2 = xr.DataArray((emis_points2_grouped['x'].values.astype(int)).tolist(),dims=["d"])
                    ind_y2 = xr.DataArray((emis_points2_grouped['y'].values.astype(int)).tolist(),dims=["d"])
                    ind_t2 = xr.DataArray((emis_points2_grouped['Time'].values.astype(int)).tolist(),dims=["d"])
                    ind_bottom_top2 = xr.DataArray((emis_points2_grouped['bottom_top'].values.astype(int)).tolist(),dims=["d"])
                    
                    emis1[variable].data[ind_t1,ind_bottom_top1,ind_y1,ind_x1] = emis_points1_grouped[variable].astype(float).values
                    emis2[variable].data[ind_t2,ind_bottom_top2,ind_y2,ind_x2] = emis_points2_grouped[variable].astype(float).values
        
        if inputs.gridded_files:
            if variable in area_emis1.data_vars:

                # add in the area emissions if the variable has any data
                if area_emis_sum_variables[variable]>0:
                    emis1[variable].loc[dict(bottom_top=0)]= emis1[variable].loc[dict(bottom_top=0)] + area_emis1[variable]*unit_conv 
                    emis2[variable].loc[dict(bottom_top=0)]= emis2[variable].loc[dict(bottom_top=0)] + area_emis2[variable]*unit_conv 
        
        # set variable attributes 
        emis1[variable].attrs = {}
        emis1[variable].attrs['units']=out_nc_file.variables[variable]['units']
        emis1[variable].attrs['description']=out_nc_file.variables[variable]['description']
        emis2[variable].attrs = {}
        emis2[variable].attrs['units']=out_nc_file.variables[variable]['units']
        emis2[variable].attrs['description']=out_nc_file.variables[variable]['description']
        
        # Set encoding
        encoding_dict = {'dtype': 'float32', 'chunksizes':(1,1,emis1.sizes['south_north'], emis1.sizes['west_east']),
                      'zlib': True, 'complevel': 1, '_FillValue': None}
        emis1[variable].encoding=encoding_dict
        emis1['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}
        emis1['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}
        #emis1['Times'].encoding={'char_dim_name':'DateStrLen'}
        emis2[variable].encoding=encoding_dict
        emis2['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}
        emis2['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}
        #emis1['Times'].encoding={'char_dim_name':'DateStrLen'}
        
        if dow == 'weekdy':
            emis1.attrs['day_of_week']='weekday'
            emis2.attrs['day_of_week']='weekday'
        elif dow == 'satdy':
            emis1.attrs['day_of_week']='saturday'
            emis2.attrs['day_of_week']='saturday'
        elif dow== 'sundy':
            emis1.attrs['day_of_week']='sunday'
            emis2.attrs['day_of_week']='sunday'
            
        emis1.attrs['Sector'] = inputs.sector_title_output
        emis2.attrs['Sector'] = inputs.sector_title_output
        
        if c ==0:
            emis1.to_netcdf(inputs.out_file1.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow),format='netCDF4',engine='netcdf4')
            emis2.to_netcdf(inputs.out_file2.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow),format='netCDF4',engine='netcdf4')
        else:
            emis1.to_netcdf(inputs.out_file1.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow),format='netCDF4',engine='netcdf4',mode='a')
            emis2.to_netcdf(inputs.out_file2.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow),format='netCDF4',engine='netcdf4',mode='a')
        
        emis1  = emis1.drop_vars([variable])
        emis2  = emis2.drop_vars([variable])
        
        if inputs.gridded_files:
            del area_emis1,area_emis2
        if inputs.point_files:
            del point_emis1,point_emis2
        
        time_loop_variable_end = time.time()
        print('Total elapsed time for variable loop (s): ' + str(round(time_loop_variable_end-time_loop_variable_start)))

def aggregate_by_xy(emis,var,sizes):
    
    e = emis.to_dataframe(name=var).reset_index()[['Time','bottom_top','x','y',var]].groupby(['Time','bottom_top','x','y']).sum().reset_index(['x','y'])
    e['id'] = e.groupby(['x','y']).ngroup()
    e.reset_index(inplace=True)#.set_index('id',append=True,inplace=True)
    
    e.drop(e.loc[e['x']>sizes['west_east']-1].index,inplace=True)
    e.drop(e.loc[e['x']<0].index,inplace=True)
    e.drop(e.loc[e['y']>sizes['south_north']-1].index,inplace=True)
    e.drop(e.loc[e['y']<0].index,inplace=True)
    
    return e
    
def main():
    #pd.options.mode.chained_assignment = 'raise'
    start_time = time.time()
    print('START OF PROGRAM: ')
    print('--------------------------------------------------')
    
    wrfchemi_ds = xr.open_dataset(inputs.wrfchemi_example,cache=False)
    # calc height ASL
    wrfchemi_ds['h_asl']= (wrfchemi_ds['PHB'] + wrfchemi_ds['PH'])/9.81 # height ASL (m)
    wrfchemi_ds['h_agl']= wrfchemi_ds['h_asl']-wrfchemi_ds['h_asl'].loc[dict(bottom_top_stag=0)]
    
    out_ds_empty = make_output_ds_structure(wrfchemi_ds)
    
    for year in inputs.years:
        for month in inputs.months:
            for dow in inputs.dows:
                
                # this needs to be done each month since number of points can change
                if inputs.point_files:
                    frac1,frac2,x_y = calc_fracs(wrfchemi_ds=wrfchemi_ds,file=inputs.point_files[0].replace('[YYYY]',str(year)).replace('[YY]',str(year)[2:]).replace('[MM]',str(month).zfill(2)).replace('[DOW]',dow).replace('[dayhalf]','00to12Z'))
                else:
                    frac1=None
                    frac2=None
                    x_y=None
                
                place_vars(year,month,dow,out_ds_empty,frac1,frac2,x_y)
    
    
    print('END OF PROGRAM: ')
    print('--------------------------------------------------')
    end_time = time.time()
    print('Total elapsed time (s): ' + str(round(end_time-start_time)))
    
# End Main

if __name__ == "__main__":
    main()

    



    

