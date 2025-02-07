#!/usr/bin/env python
# coding: utf-8

# In[1]:


#mamba activate /home/clyu/miniconda3/envs/xesmf_env
from netCDF4 import Dataset
import numpy as np
import reverse_geocoder as rg
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
mpl.use('Agg')
#from mpl_toolkits.basemap import Basemap
import matplotlib.ticker
from matplotlib.patches import Polygon
import xarray as xr
#import xesmf as xe
import pandas as pd
#import ESMF
import pyproj
import os
import scipy.stats as stats
import pyreadr
import reverse_geocoder as rg
from rpy2.robjects.packages import importr
base = importr('base')
import statistics


# In[2]:


version = 'V7_GRA2PES2022'
mm = '12'
mm_index = 12
#base_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_wrfchemi2022_input/point/Month'+mm+'/TotlPoint_OG'
#append_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_wrfchemi2022_input/point_append_extra/Month'+mm+'/TotlPoint_OG'


base_dir = '/wrk/users/charkins/emissions/V7_GRA2PES/POINT22_ncf/Month'+mm+'/PtOG'
append_dir = '/wrk/users/charkins/emissions/V7_GRA2PES/POINT22_ncf_append_extra/Month'+mm+'/PtOG'

# In[3]:


#read GHGRP OG points
GHGRP_ng_proc_all_stackinfo = pd.read_csv("/wrk/csd4/clyu/GHG_CO2/Improving_inventory/"+version+"/extra_points/input/GHGRP_ng_proc_2022_FCPE_all_stackinfo_emis.csv")
LON_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,1]
LAT_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,2]
STATE_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,3]
CO2_stack_FC_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,4] #metric tons/year
CO2_stack_PE_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,5] #metric tons/year
CH4_stack_FC_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,6] #metric tons/year
CH4_stack_PE_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,7] #metric tons/year
ERPTYPE_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,8]
STKHGT_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,9]
STKDIAM_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,10]
STKTEMP_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,11]
STKFLOW_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,12]
STKVEL_ng_proc = GHGRP_ng_proc_all_stackinfo.iloc[:,13]


# In[4]:


#Calculate state-level AQ species to ffCO2 emission ratios


# In[5]:


#OG
#combine FC PE
CO2_stack_FCPE_ng_proc = CO2_stack_FC_ng_proc + CO2_stack_PE_ng_proc
CH4_stack_FCPE_ng_proc = CH4_stack_FC_ng_proc + CH4_stack_PE_ng_proc

#get CH4 to CO2 FCPE ratio from GHGRP, only available for selected states that has GHGRP points
states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

process_vector = ['OnG']

###########################################################################################################
#ng_proc
#sum CO2 and CH4 respectively at state-level
CO2_stack_FCPE_ng_proc_STs_sum = np.zeros([len(states_abb_vector)])
CH4_stack_FCPE_ng_proc_STs_sum = np.zeros([len(states_abb_vector)])

STATE_ng_proc_uniq = np.unique(STATE_ng_proc)
STATE_ng_proc = np.array(STATE_ng_proc)

staten = 0
for ST in states_abb_vector:
    if ST in STATE_ng_proc_uniq:
        ST_ind, = np.where(STATE_ng_proc==ST)
        CO2_stack_FCPE_ng_proc_ST = CO2_stack_FCPE_ng_proc[ST_ind]
        CO2_stack_FCPE_ng_proc_ST_sum = np.nansum(CO2_stack_FCPE_ng_proc_ST)
        CO2_stack_FCPE_ng_proc_STs_sum[staten] = CO2_stack_FCPE_ng_proc_ST_sum
        CH4_stack_FCPE_ng_proc_ST = CH4_stack_FCPE_ng_proc[ST_ind]
        CH4_stack_FCPE_ng_proc_ST_sum = np.nansum(CH4_stack_FCPE_ng_proc_ST)
        CH4_stack_FCPE_ng_proc_STs_sum[staten] = CH4_stack_FCPE_ng_proc_ST_sum
    staten += 1

CO2_stack_FCPE_ng_proc_national = np.nansum(CO2_stack_FCPE_ng_proc_STs_sum)
CH4_stack_FCPE_ng_proc_national = np.nansum(CH4_stack_FCPE_ng_proc_STs_sum)

#if state-level summed CO2 or CH4 is zero, patch both CO2 and CH4 with national-level summed CO2 and CH4
staten = 0
for ST in states_abb_vector:
    if CO2_stack_FCPE_ng_proc_STs_sum[staten] == 0 or CH4_stack_FCPE_ng_proc_STs_sum[staten] == 0:
        CO2_stack_FCPE_ng_proc_STs_sum[staten] = CO2_stack_FCPE_ng_proc_national
        CH4_stack_FCPE_ng_proc_STs_sum[staten] = CH4_stack_FCPE_ng_proc_national
    staten += 1

#calculate state-level CH4/CO2 ratio using state-level sum data with national-level sum patch
ratio_CH4_CO2_FCPE_ng_proc_STs = CH4_stack_FCPE_ng_proc_STs_sum/CO2_stack_FCPE_ng_proc_STs_sum
print("ratio_CH4_CO2_FCPE_ng_proc_STs",ratio_CH4_CO2_FCPE_ng_proc_STs)

###########################################################################################################
proc_state_emisCH4dCO2_FCPE = np.empty([len(process_vector),len(states_vector)], dtype=object)
proc_state_emisCH4dCO2_FCPE[0,:] = ratio_CH4_CO2_FCPE_ng_proc_STs

###########################################################################################################
#get process-specific other AQ to CH4 ratio from
#/wrk/csd4/clyu/GHG_CO2/Improving_inventory/'+version+'/scal_sum_ncf/base2017_rds/point/Month00/PtOnG/weekdy

proc_spec_state_emis_OG = np.empty([len(process_vector),len(species_vector),len(states_vector)], dtype=object)

procn=0
for proc in process_vector:
    
    specn = 0
    for spec in species_vector:
        PtOG_spec_weekdy = pyreadr.read_r('/wrk/csd4/clyu/GHG_CO2/Improving_inventory/'+version+'/scal_sum_ncf/base2017_rds/point/Month00/Pt'+proc+'/weekdy/Pt'+proc+'_'+spec+'_weekdy.rds')
        PtOG_spec_weekdy_np = np.array(PtOG_spec_weekdy[None])
        STATE_col = PtOG_spec_weekdy_np[:,2]
        dayav_col = PtOG_spec_weekdy_np[:,4] #metric tons per day

        staten = 0
        for state in states_vector:
            dayav_col_state = dayav_col[np.where(STATE_col==state)]
            dayav_col_state_sum = np.nansum(dayav_col_state)
            proc_spec_state_emis_OG[procn,specn,staten] = dayav_col_state_sum
            staten += 1
        specn += 1
    procn += 1

print("proc_spec_state_emis_OG",proc_spec_state_emis_OG)

#if CH4 is 0, make it nan
for pp in range(0,procn):
    for st in range(0,staten):
        if proc_spec_state_emis_OG[pp,7,st] == 0:
            proc_spec_state_emis_OG[pp,7,st] = np.nan

proc_spec_state_emisXdCH4_OG = np.empty([len(process_vector),len(species_vector),len(states_vector)], dtype=object)

for pp in range(0,procn):
    for sp in range(0,specn):
        for st in range(0,staten):
            proc_spec_state_emisXdCH4_OG[pp,sp,st] = proc_spec_state_emis_OG[pp,sp,st]/proc_spec_state_emis_OG[pp,7,st]

print("proc_spec_state_emisXdCH4_OG",proc_spec_state_emisXdCH4_OG)

proc_spec_state_emisXdCO2_OG = np.empty([len(process_vector),len(species_vector),len(states_vector)])

for pp in range(0,procn):
    for sp in range(0,specn):
        for st in range(0,staten):
            proc_spec_state_emisXdCO2_OG[pp,sp,st] = proc_spec_state_emisXdCH4_OG[pp,sp,st] * proc_state_emisCH4dCO2_FCPE[pp,st]

#make nan to zero
proc_spec_state_emisXdCO2_OG[np.isnan(proc_spec_state_emisXdCO2_OG)] = 0.0
print("proc_spec_state_emisXdCO2_OG",proc_spec_state_emisXdCO2_OG)


# In[6]:


#Mimic Stu's RELPT program, convert 2022 annual total emissions 
#to day of week averaged daily total, and 24 hours emissions


# In[7]:


#First, getting d.o.w fractions from RELPT output emission files


# In[8]:


#OG, process-specific
#read processed emissions that has weekdy, satdy, sundy to get d.o.w. fractions
###########################################################################################################################
PtOnG_CO_weekdy = pyreadr.read_r('/wrk/csd4/clyu/GHG_CO2/Improving_inventory/'+version+'/scal_sum_ncf/base2017_rds/point/Month00/PtOnG/weekdy/PtOnG_CO_weekdy.rds')
PtOnG_CO_satdy = pyreadr.read_r('/wrk/csd4/clyu/GHG_CO2/Improving_inventory/'+version+'/scal_sum_ncf/base2017_rds/point/Month00/PtOnG/satdy/PtOnG_CO_satdy.rds')
PtOnG_CO_sundy = pyreadr.read_r('/wrk/csd4/clyu/GHG_CO2/Improving_inventory/'+version+'/scal_sum_ncf/base2017_rds/point/Month00/PtOnG/sundy/PtOnG_CO_sundy.rds')

PtOnG_CO_weekdy_np = np.array(PtOnG_CO_weekdy[None])
PtOnG_CO_weekdy_dayav = PtOnG_CO_weekdy_np[:,4] #metric tons per day
PtOnG_CO_weekdy_dayav_sum = np.nansum(PtOnG_CO_weekdy_dayav) #metric tons per day in the whole domain
print("PtOnG_CO_weekdy_dayav_sum",PtOnG_CO_weekdy_dayav_sum)

PtOnG_CO_satdy_np = np.array(PtOnG_CO_satdy[None])
PtOnG_CO_satdy_dayav = PtOnG_CO_satdy_np[:,4] #metric tons per day
PtOnG_CO_satdy_dayav_sum = np.nansum(PtOnG_CO_satdy_dayav) #metric tons per day in the whole domain
print("PtOnG_CO_satdy_dayav_sum",PtOnG_CO_satdy_dayav_sum)

PtOnG_CO_sundy_np = np.array(PtOnG_CO_sundy[None])
PtOnG_CO_sundy_dayav = PtOnG_CO_sundy_np[:,4] #metric tons per day
PtOnG_CO_sundy_dayav_sum = np.nansum(PtOnG_CO_sundy_dayav) #metric tons per day in the whole domain
print("PtOnG_CO_sundy_dayav_sum",PtOnG_CO_sundy_dayav_sum)

#metric tons per week in the whole domain
PtOnG_CO_week_dayav_sum = PtOnG_CO_weekdy_dayav_sum*5 + PtOnG_CO_satdy_dayav_sum + PtOnG_CO_sundy_dayav_sum
PtOnG_CO_5weekdy_fraction = PtOnG_CO_weekdy_dayav_sum*5/PtOnG_CO_week_dayav_sum
PtOnG_CO_1satdy_fraction = PtOnG_CO_satdy_dayav_sum/PtOnG_CO_week_dayav_sum
PtOnG_CO_1sundy_fraction = PtOnG_CO_sundy_dayav_sum/PtOnG_CO_week_dayav_sum
print("PtOnG_CO_5weekdy_fraction",PtOnG_CO_5weekdy_fraction)
print("PtOnG_CO_1satdy_fraction",PtOnG_CO_1satdy_fraction)
print("PtOnG_CO_1sundy_fraction",PtOnG_CO_1sundy_fraction)


# In[9]:


#Apply fuel or process-specific d.o.w. fractions to annual total emissions


# In[10]:


#OG, apply process-specific d.o.w. fractions
num_weekdys = 5*52
num_satdys = 52
num_sundys = 52

###################################################################################################
#ng_proc
###################################################################################################
#CO2
###################################################################################################

#FCPE
dayav_CO2_FCPE_MetricTon_2022_ng_proc_weekdy = np.zeros(len(LON_ng_proc))
dayav_CO2_FCPE_MetricTon_2022_ng_proc_satdy = np.zeros(len(LON_ng_proc))
dayav_CO2_FCPE_MetricTon_2022_ng_proc_sundy = np.zeros(len(LON_ng_proc))

for pt in range(0,len(LON_ng_proc)):
    dayav_CO2_FCPE_MetricTon_2022_ng_proc_weekdy[pt] = CO2_stack_FCPE_ng_proc[pt]*PtOnG_CO_5weekdy_fraction/num_weekdys
    dayav_CO2_FCPE_MetricTon_2022_ng_proc_satdy[pt] = CO2_stack_FCPE_ng_proc[pt]*PtOnG_CO_1satdy_fraction/num_satdys
    dayav_CO2_FCPE_MetricTon_2022_ng_proc_sundy[pt] = CO2_stack_FCPE_ng_proc[pt]*PtOnG_CO_1sundy_fraction/num_sundys

###################################################################################################
#CH4
###################################################################################################

#FCPE
dayav_CH4_FCPE_MetricTon_2022_ng_proc_weekdy = np.zeros(len(LON_ng_proc))
dayav_CH4_FCPE_MetricTon_2022_ng_proc_satdy = np.zeros(len(LON_ng_proc))
dayav_CH4_FCPE_MetricTon_2022_ng_proc_sundy = np.zeros(len(LON_ng_proc))

for pt in range(0,len(LON_ng_proc)):
    dayav_CH4_FCPE_MetricTon_2022_ng_proc_weekdy[pt] = CH4_stack_FCPE_ng_proc[pt]*PtOnG_CO_5weekdy_fraction/num_weekdys
    dayav_CH4_FCPE_MetricTon_2022_ng_proc_satdy[pt] = CH4_stack_FCPE_ng_proc[pt]*PtOnG_CO_1satdy_fraction/num_satdys
    dayav_CH4_FCPE_MetricTon_2022_ng_proc_sundy[pt] = CH4_stack_FCPE_ng_proc[pt]*PtOnG_CO_1sundy_fraction/num_sundys


# In[11]:


#get fuel or process and d.o.w.-specific state-level 24 hours temporal profile


# In[12]:


#OG, process and d.o.w.-specific state-level 24 hours temporal profile
###################################################################################################
processes_vector = ['OnG']

dow_vector = ['weekdy','satdy','sundy']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

ETstates = ['Maine','New Hampshire','Vermont','Massachusetts','Rhode Island','Connecticut','New York','New Jersey',
            'Pennsylvania','Delaware','Maryland','District of Columbia','Virginia','West Virginia','North Carolina',
            'South Carolina','Georgia','Florida','Ohio','Michigan', 'Indiana','Kentucky']
CTstates = ['Alabama','Arkansas','Illinois','Iowa','Kansas','Louisiana','Minnesota',
            'Mississippi','Missouri','Nebraska','North Dakota','Oklahoma','South Dakota',
            'Texas','Tennessee','Wisconsin']
MTstates = ['Arizona','Colorado','Idaho','Montana','New Mexico','Utah','Wyoming']
PTstates = ['California','Washington','Oregon','Nevada']

for proc in processes_vector:
    for dow in dow_vector:
        print("proc",proc)
        print("dow",dow)
        PtOG_proc_dow = pyreadr.read_r('/wrk/csd4/clyu/GHG_CO2/Improving_inventory/'+version+'/scal_sum_ncf/base2017_rds/point/Month00/Pt'+proc+'/'+dow+'/Pt'+proc+'_CO_'+dow+'.rds')
        PtOG_proc_dow_np = np.array(PtOG_proc_dow[None])
        STATE_col = PtOG_proc_dow_np[:,2]
        PtOG_proc_dow_dayav = PtOG_proc_dow_np[:,4] #metric tons per day
        PtOG_proc_dow_states_HRall_frac = np.zeros([len(states_vector),24])

        #prepare time zone level 24 hr profile for states have no profile
        #ET
        PtOG_proc_dow_ETstates_HRall_frac = np.zeros([24])
        PtOG_proc_dow_dayav_ETstate = PtOG_proc_dow_dayav[np.where(np.isin(STATE_col, ETstates))]
        PtOG_proc_dow_dayav_ETstate_sum = np.nansum(PtOG_proc_dow_dayav_ETstate)
        
        for hh in range(0,24):
            PtOG_proc_dow_hh = PtOG_proc_dow_np[:,hh+5] #metric tons in this hour
            PtOG_proc_dow_hh_ETstate = PtOG_proc_dow_hh[np.where(np.isin(STATE_col, ETstates))]
            PtOG_proc_dow_hh_ETstate_sum = np.nansum(PtOG_proc_dow_hh_ETstate)
            PtOG_proc_dow_ETstates_HRall_frac[hh] = PtOG_proc_dow_hh_ETstate_sum/PtOG_proc_dow_dayav_ETstate_sum
        print("PtOG_proc_dow_ETstates_HRall_frac",PtOG_proc_dow_ETstates_HRall_frac)

        #CT
        PtOG_proc_dow_CTstates_HRall_frac = np.zeros([24])
        PtOG_proc_dow_dayav_CTstate = PtOG_proc_dow_dayav[np.where(np.isin(STATE_col, CTstates))]
        PtOG_proc_dow_dayav_CTstate_sum = np.nansum(PtOG_proc_dow_dayav_CTstate)
        
        for hh in range(0,24):
            PtOG_proc_dow_hh = PtOG_proc_dow_np[:,hh+5] #mCTric tons in this hour
            PtOG_proc_dow_hh_CTstate = PtOG_proc_dow_hh[np.where(np.isin(STATE_col, CTstates))]
            PtOG_proc_dow_hh_CTstate_sum = np.nansum(PtOG_proc_dow_hh_CTstate)
            PtOG_proc_dow_CTstates_HRall_frac[hh] = PtOG_proc_dow_hh_CTstate_sum/PtOG_proc_dow_dayav_CTstate_sum
        print("PtOG_proc_dow_CTstates_HRall_frac",PtOG_proc_dow_CTstates_HRall_frac)

        #MT
        PtOG_proc_dow_MTstates_HRall_frac = np.zeros([24])
        PtOG_proc_dow_dayav_MTstate = PtOG_proc_dow_dayav[np.where(np.isin(STATE_col, MTstates))]
        PtOG_proc_dow_dayav_MTstate_sum = np.nansum(PtOG_proc_dow_dayav_MTstate)
        
        for hh in range(0,24):
            PtOG_proc_dow_hh = PtOG_proc_dow_np[:,hh+5] #mMTric tons in this hour
            PtOG_proc_dow_hh_MTstate = PtOG_proc_dow_hh[np.where(np.isin(STATE_col, MTstates))]
            PtOG_proc_dow_hh_MTstate_sum = np.nansum(PtOG_proc_dow_hh_MTstate)
            PtOG_proc_dow_MTstates_HRall_frac[hh] = PtOG_proc_dow_hh_MTstate_sum/PtOG_proc_dow_dayav_MTstate_sum
        print("PtOG_proc_dow_MTstates_HRall_frac",PtOG_proc_dow_MTstates_HRall_frac)

        #PT
        PtOG_proc_dow_PTstates_HRall_frac = np.zeros([24])
        PtOG_proc_dow_dayav_PTstate = PtOG_proc_dow_dayav[np.where(np.isin(STATE_col, PTstates))]
        PtOG_proc_dow_dayav_PTstate_sum = np.nansum(PtOG_proc_dow_dayav_PTstate)
        
        for hh in range(0,24):
            PtOG_proc_dow_hh = PtOG_proc_dow_np[:,hh+5] #mPTric tons in this hour
            PtOG_proc_dow_hh_PTstate = PtOG_proc_dow_hh[np.where(np.isin(STATE_col, PTstates))]
            PtOG_proc_dow_hh_PTstate_sum = np.nansum(PtOG_proc_dow_hh_PTstate)
            PtOG_proc_dow_PTstates_HRall_frac[hh] = PtOG_proc_dow_hh_PTstate_sum/PtOG_proc_dow_dayav_PTstate_sum
        print("PtOG_proc_dow_PTstates_HRall_frac",PtOG_proc_dow_PTstates_HRall_frac)
        
        #Get state-level 24 hrs profile
        staten = 0
        for state in states_vector:
            print("state",state)
            PtOG_proc_dow_dayav_state = PtOG_proc_dow_dayav[np.where(STATE_col==state)]
            PtOG_proc_dow_dayav_state_sum = np.nansum(PtOG_proc_dow_dayav_state)
            
            if PtOG_proc_dow_dayav_state_sum > 0:
            
                for hh in range(0,24):
                    PtOG_proc_dow_hh = PtOG_proc_dow_np[:,hh+5] #metric tons in this hour

                    PtOG_proc_dow_hh_state = PtOG_proc_dow_hh[np.where(STATE_col==state)]
                    PtOG_proc_dow_hh_state_sum = np.nansum(PtOG_proc_dow_hh_state)
                    PtOG_proc_dow_states_HRall_frac[staten,hh] = PtOG_proc_dow_hh_state_sum/PtOG_proc_dow_dayav_state_sum
            
            else:
                if state in ETstates:
                    PtOG_proc_dow_states_HRall_frac[staten,:] = PtOG_proc_dow_ETstates_HRall_frac
                elif state in CTstates:
                    PtOG_proc_dow_states_HRall_frac[staten,:] = PtOG_proc_dow_CTstates_HRall_frac
                elif state in MTstates:
                    PtOG_proc_dow_states_HRall_frac[staten,:] = PtOG_proc_dow_MTstates_HRall_frac
                elif state in PTstates:
                    PtOG_proc_dow_states_HRall_frac[staten,:] = PtOG_proc_dow_PTstates_HRall_frac

            #sanity check
            PtOG_proc_dow_state_HRall_frac = np.nansum(PtOG_proc_dow_states_HRall_frac[staten,:])
            print("PtOG_proc_dow_state_HRall_frac",PtOG_proc_dow_state_HRall_frac)
            
            staten += 1
                    
        if proc == 'OnG':
            if dow == 'weekdy':
                PtOnG_weekdy_states_HRall_frac = PtOG_proc_dow_states_HRall_frac
            elif dow == 'satdy':
                PtOnG_satdy_states_HRall_frac = PtOG_proc_dow_states_HRall_frac
            elif dow == 'sundy':
                PtOnG_sundy_states_HRall_frac = PtOG_proc_dow_states_HRall_frac


# In[13]:


#Apply 24 hours profile


# In[14]:


#OG, process and d.o.w.-specific state-level 24 hours temporal profile
###########################################################################################################################
states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#ng_proc
###################################################################################################
#CO2
###################################################################################################
HRall_CO2_FCPE_MetricTon_2022_ng_proc_weekdy = np.zeros([24,len(LON_ng_proc)])
HRall_CO2_FCPE_MetricTon_2022_ng_proc_satdy = np.zeros([24,len(LON_ng_proc)])
HRall_CO2_FCPE_MetricTon_2022_ng_proc_sundy = np.zeros([24,len(LON_ng_proc)])

for pt in range(0,len(LON_ng_proc)):
    #print("pt",pt)
    state_cur = STATE_ng_proc[pt]
    if state_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_cur)
    #print("state_cur",state_cur)
    #print("state_index",state_index)
    
    HRall_CO2_FCPE_MetricTon_2022_ng_proc_weekdy[:,pt]= dayav_CO2_FCPE_MetricTon_2022_ng_proc_weekdy[pt]*PtOnG_weekdy_states_HRall_frac[state_index,:]
    HRall_CO2_FCPE_MetricTon_2022_ng_proc_satdy[:,pt]= dayav_CO2_FCPE_MetricTon_2022_ng_proc_satdy[pt]*PtOnG_satdy_states_HRall_frac[state_index,:]
    HRall_CO2_FCPE_MetricTon_2022_ng_proc_sundy[:,pt]= dayav_CO2_FCPE_MetricTon_2022_ng_proc_sundy[pt]*PtOnG_sundy_states_HRall_frac[state_index,:]

print("dayav_CO2_FCPE_MetricTon_2022_ng_proc_weekdy",np.nansum(dayav_CO2_FCPE_MetricTon_2022_ng_proc_weekdy))
print("HRall_CO2_FCPE_MetricTon_2022_ng_proc_weekdy",np.nansum(HRall_CO2_FCPE_MetricTon_2022_ng_proc_weekdy))
print("dayav_CO2_FCPE_MetricTon_2022_ng_proc_satdy",np.nansum(dayav_CO2_FCPE_MetricTon_2022_ng_proc_satdy))
print("HRall_CO2_FCPE_MetricTon_2022_ng_proc_satdy",np.nansum(HRall_CO2_FCPE_MetricTon_2022_ng_proc_satdy))
print("dayav_CO2_FCPE_MetricTon_2022_ng_proc_sundy",np.nansum(dayav_CO2_FCPE_MetricTon_2022_ng_proc_sundy))
print("HRall_CO2_FCPE_MetricTon_2022_ng_proc_sundy",np.nansum(HRall_CO2_FCPE_MetricTon_2022_ng_proc_sundy))

###################################################################################################
#CH4
###################################################################################################
HRall_CH4_FCPE_MetricTon_2022_ng_proc_weekdy = np.zeros([24,len(LON_ng_proc)])
HRall_CH4_FCPE_MetricTon_2022_ng_proc_satdy = np.zeros([24,len(LON_ng_proc)])
HRall_CH4_FCPE_MetricTon_2022_ng_proc_sundy = np.zeros([24,len(LON_ng_proc)])

for pt in range(0,len(LON_ng_proc)):
    #print("pt",pt)
    state_cur = STATE_ng_proc[pt]
    if state_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_cur)
    #print("state_cur",state_cur)
    #print("state_index",state_index)
    
    HRall_CH4_FCPE_MetricTon_2022_ng_proc_weekdy[:,pt]= dayav_CH4_FCPE_MetricTon_2022_ng_proc_weekdy[pt]*PtOnG_weekdy_states_HRall_frac[state_index,:]
    HRall_CH4_FCPE_MetricTon_2022_ng_proc_satdy[:,pt]= dayav_CH4_FCPE_MetricTon_2022_ng_proc_satdy[pt]*PtOnG_satdy_states_HRall_frac[state_index,:]
    HRall_CH4_FCPE_MetricTon_2022_ng_proc_sundy[:,pt]= dayav_CH4_FCPE_MetricTon_2022_ng_proc_sundy[pt]*PtOnG_sundy_states_HRall_frac[state_index,:]

print("dayav_CH4_FCPE_MetricTon_2022_ng_proc_weekdy",np.nansum(dayav_CH4_FCPE_MetricTon_2022_ng_proc_weekdy))
print("HRall_CH4_FCPE_MetricTon_2022_ng_proc_weekdy",np.nansum(HRall_CH4_FCPE_MetricTon_2022_ng_proc_weekdy))
print("dayav_CH4_FCPE_MetricTon_2022_ng_proc_satdy",np.nansum(dayav_CH4_FCPE_MetricTon_2022_ng_proc_satdy))
print("HRall_CH4_FCPE_MetricTon_2022_ng_proc_satdy",np.nansum(HRall_CH4_FCPE_MetricTon_2022_ng_proc_satdy))
print("dayav_CH4_FCPE_MetricTon_2022_ng_proc_sundy",np.nansum(dayav_CH4_FCPE_MetricTon_2022_ng_proc_sundy))
print("HRall_CH4_FCPE_MetricTon_2022_ng_proc_sundy",np.nansum(HRall_CH4_FCPE_MetricTon_2022_ng_proc_sundy))


# In[15]:


#scale ng_proc from 2022 annual average in metric tons/hr to 2022mm
PtOG_monthly = pd.read_csv("/wrk/csd4/charkins/emissions/GRA2PES/V7_NRT_scaling/POINT22_202411/input/PtOG_monthly.csv")

PtOnG_sf = PtOG_monthly.iloc[:,1]

#get 2022 four season scaling from 2017 annual average
PtOnG_sf_jan = PtOnG_sf[(1-1)]
PtOnG_sf_apr = PtOnG_sf[(4-1)]
PtOnG_sf_jul = PtOnG_sf[(7-1)]
PtOnG_sf_oct = PtOnG_sf[(10-1)]

#calculate 2022 annual average scaling from 2017 annual average
PtOnG_sf_aavg = (PtOnG_sf_jan+PtOnG_sf_apr+PtOnG_sf_jul+PtOnG_sf_oct)/4

#get 2022mm scaling from 2017 annual average
PtOnG_sf_mm = PtOnG_sf[(mm_index-1)]

#calculate 2022mm scaling from 2022 annual average
PtOnG_sf_2022mm = PtOnG_sf_mm/PtOnG_sf_aavg

print("PtOnG_sf_2022mm",PtOnG_sf_2022mm)

HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_weekdy = HRall_CO2_FCPE_MetricTon_2022_ng_proc_weekdy * PtOnG_sf_2022mm
HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_satdy = HRall_CO2_FCPE_MetricTon_2022_ng_proc_satdy * PtOnG_sf_2022mm
HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_sundy = HRall_CO2_FCPE_MetricTon_2022_ng_proc_sundy * PtOnG_sf_2022mm

HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_weekdy = HRall_CH4_FCPE_MetricTon_2022_ng_proc_weekdy * PtOnG_sf_2022mm
HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_satdy = HRall_CH4_FCPE_MetricTon_2022_ng_proc_satdy * PtOnG_sf_2022mm
HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_sundy = HRall_CH4_FCPE_MetricTon_2022_ng_proc_sundy * PtOnG_sf_2022mm


# In[16]:


#append extra points to original point file that is input to wrfchemi assembly program

###################################################################################################
#weekdy, 00to12Z

###################################################################################################
#read original variables
TotlPoint_00to12Z_weekdy_fn = base_dir+'/weekdy/PtOG_00to12Z.nc'
TotlPoint_00to12Z_weekdy_file = Dataset(TotlPoint_00to12Z_weekdy_fn,mode='r',open=True)
TotlPoint_00to12Z_weekdy_ITYPE = TotlPoint_00to12Z_weekdy_file.variables['ITYPE'][:]
TotlPoint_00to12Z_weekdy_STKht = TotlPoint_00to12Z_weekdy_file.variables['STKht'][:]
TotlPoint_00to12Z_weekdy_STKdiam = TotlPoint_00to12Z_weekdy_file.variables['STKdiam'][:]
TotlPoint_00to12Z_weekdy_STKtemp = TotlPoint_00to12Z_weekdy_file.variables['STKtemp'][:]
TotlPoint_00to12Z_weekdy_STKve = TotlPoint_00to12Z_weekdy_file.variables['STKve'][:]
TotlPoint_00to12Z_weekdy_STKflw = TotlPoint_00to12Z_weekdy_file.variables['STKflw'][:]
TotlPoint_00to12Z_weekdy_FUGht = TotlPoint_00to12Z_weekdy_file.variables['FUGht'][:]
TotlPoint_00to12Z_weekdy_XLONG = TotlPoint_00to12Z_weekdy_file.variables['XLONG'][:]
TotlPoint_00to12Z_weekdy_XLAT = TotlPoint_00to12Z_weekdy_file.variables['XLAT'][:]
TotlPoint_00to12Z_weekdy_CO2 = TotlPoint_00to12Z_weekdy_file.variables['CO2'][:][:] 
TotlPoint_00to12Z_weekdy_CO = TotlPoint_00to12Z_weekdy_file.variables['CO'][:][:] 
TotlPoint_00to12Z_weekdy_NH3 = TotlPoint_00to12Z_weekdy_file.variables['NH3'][:][:] 
TotlPoint_00to12Z_weekdy_NOX = TotlPoint_00to12Z_weekdy_file.variables['NOX'][:][:] 
TotlPoint_00to12Z_weekdy_PM10_PRI = TotlPoint_00to12Z_weekdy_file.variables['PM10-PRI'][:][:] 
TotlPoint_00to12Z_weekdy_PM25_PRI = TotlPoint_00to12Z_weekdy_file.variables['PM25-PRI'][:][:] 
TotlPoint_00to12Z_weekdy_SO2 = TotlPoint_00to12Z_weekdy_file.variables['SO2'][:][:] 
TotlPoint_00to12Z_weekdy_VOC = TotlPoint_00to12Z_weekdy_file.variables['VOC'][:][:] 
TotlPoint_00to12Z_weekdy_HC01 = TotlPoint_00to12Z_weekdy_file.variables['HC01'][:][:] 
TotlPoint_00to12Z_weekdy_HC02 = TotlPoint_00to12Z_weekdy_file.variables['HC02'][:][:] 
TotlPoint_00to12Z_weekdy_HC03 = TotlPoint_00to12Z_weekdy_file.variables['HC03'][:][:] 
TotlPoint_00to12Z_weekdy_HC04 = TotlPoint_00to12Z_weekdy_file.variables['HC04'][:][:] 
TotlPoint_00to12Z_weekdy_HC05 = TotlPoint_00to12Z_weekdy_file.variables['HC05'][:][:] 
TotlPoint_00to12Z_weekdy_HC06 = TotlPoint_00to12Z_weekdy_file.variables['HC06'][:][:] 
TotlPoint_00to12Z_weekdy_HC07 = TotlPoint_00to12Z_weekdy_file.variables['HC07'][:][:] 
TotlPoint_00to12Z_weekdy_HC08 = TotlPoint_00to12Z_weekdy_file.variables['HC08'][:][:] 
TotlPoint_00to12Z_weekdy_HC09 = TotlPoint_00to12Z_weekdy_file.variables['HC09'][:][:] 
TotlPoint_00to12Z_weekdy_HC10 = TotlPoint_00to12Z_weekdy_file.variables['HC10'][:][:] 
TotlPoint_00to12Z_weekdy_HC11 = TotlPoint_00to12Z_weekdy_file.variables['HC11'][:][:] 
TotlPoint_00to12Z_weekdy_HC12 = TotlPoint_00to12Z_weekdy_file.variables['HC12'][:][:] 
TotlPoint_00to12Z_weekdy_HC13 = TotlPoint_00to12Z_weekdy_file.variables['HC13'][:][:] 
TotlPoint_00to12Z_weekdy_HC14 = TotlPoint_00to12Z_weekdy_file.variables['HC14'][:][:] 
TotlPoint_00to12Z_weekdy_HC15 = TotlPoint_00to12Z_weekdy_file.variables['HC15'][:][:] 
TotlPoint_00to12Z_weekdy_HC16 = TotlPoint_00to12Z_weekdy_file.variables['HC16'][:][:] 
TotlPoint_00to12Z_weekdy_HC17 = TotlPoint_00to12Z_weekdy_file.variables['HC17'][:][:] 
TotlPoint_00to12Z_weekdy_HC18 = TotlPoint_00to12Z_weekdy_file.variables['HC18'][:][:] 
TotlPoint_00to12Z_weekdy_HC19 = TotlPoint_00to12Z_weekdy_file.variables['HC19'][:][:] 
TotlPoint_00to12Z_weekdy_HC20 = TotlPoint_00to12Z_weekdy_file.variables['HC20'][:][:] 
TotlPoint_00to12Z_weekdy_HC21 = TotlPoint_00to12Z_weekdy_file.variables['HC21'][:][:] 
TotlPoint_00to12Z_weekdy_HC22 = TotlPoint_00to12Z_weekdy_file.variables['HC22'][:][:] 
TotlPoint_00to12Z_weekdy_HC23 = TotlPoint_00to12Z_weekdy_file.variables['HC23'][:][:] 
TotlPoint_00to12Z_weekdy_HC24 = TotlPoint_00to12Z_weekdy_file.variables['HC24'][:][:] 
TotlPoint_00to12Z_weekdy_HC25 = TotlPoint_00to12Z_weekdy_file.variables['HC25'][:][:] 
TotlPoint_00to12Z_weekdy_HC26 = TotlPoint_00to12Z_weekdy_file.variables['HC26'][:][:] 
TotlPoint_00to12Z_weekdy_HC27 = TotlPoint_00to12Z_weekdy_file.variables['HC27'][:][:] 
TotlPoint_00to12Z_weekdy_HC28 = TotlPoint_00to12Z_weekdy_file.variables['HC28'][:][:] 
TotlPoint_00to12Z_weekdy_HC29 = TotlPoint_00to12Z_weekdy_file.variables['HC29'][:][:] 
TotlPoint_00to12Z_weekdy_HC30 = TotlPoint_00to12Z_weekdy_file.variables['HC30'][:][:] 
TotlPoint_00to12Z_weekdy_HC31 = TotlPoint_00to12Z_weekdy_file.variables['HC31'][:][:] 
TotlPoint_00to12Z_weekdy_HC32 = TotlPoint_00to12Z_weekdy_file.variables['HC32'][:][:] 
TotlPoint_00to12Z_weekdy_HC33 = TotlPoint_00to12Z_weekdy_file.variables['HC33'][:][:] 
TotlPoint_00to12Z_weekdy_HC34 = TotlPoint_00to12Z_weekdy_file.variables['HC34'][:][:] 
TotlPoint_00to12Z_weekdy_HC35 = TotlPoint_00to12Z_weekdy_file.variables['HC35'][:][:] 
TotlPoint_00to12Z_weekdy_HC36 = TotlPoint_00to12Z_weekdy_file.variables['HC36'][:][:] 
TotlPoint_00to12Z_weekdy_HC37 = TotlPoint_00to12Z_weekdy_file.variables['HC37'][:][:] 
TotlPoint_00to12Z_weekdy_HC38 = TotlPoint_00to12Z_weekdy_file.variables['HC38'][:][:] 
TotlPoint_00to12Z_weekdy_HC39 = TotlPoint_00to12Z_weekdy_file.variables['HC39'][:][:] 
TotlPoint_00to12Z_weekdy_HC40 = TotlPoint_00to12Z_weekdy_file.variables['HC40'][:][:] 
TotlPoint_00to12Z_weekdy_HC41 = TotlPoint_00to12Z_weekdy_file.variables['HC41'][:][:] 
TotlPoint_00to12Z_weekdy_HC42 = TotlPoint_00to12Z_weekdy_file.variables['HC42'][:][:] 
TotlPoint_00to12Z_weekdy_HC43 = TotlPoint_00to12Z_weekdy_file.variables['HC43'][:][:] 
TotlPoint_00to12Z_weekdy_HC44 = TotlPoint_00to12Z_weekdy_file.variables['HC44'][:][:] 
TotlPoint_00to12Z_weekdy_HC45 = TotlPoint_00to12Z_weekdy_file.variables['HC45'][:][:] 
TotlPoint_00to12Z_weekdy_HC46 = TotlPoint_00to12Z_weekdy_file.variables['HC46'][:][:] 
TotlPoint_00to12Z_weekdy_HC47 = TotlPoint_00to12Z_weekdy_file.variables['HC47'][:][:] 
TotlPoint_00to12Z_weekdy_HC48 = TotlPoint_00to12Z_weekdy_file.variables['HC48'][:][:] 
TotlPoint_00to12Z_weekdy_HC49 = TotlPoint_00to12Z_weekdy_file.variables['HC49'][:][:] 
TotlPoint_00to12Z_weekdy_HC50 = TotlPoint_00to12Z_weekdy_file.variables['HC50'][:][:] 
TotlPoint_00to12Z_weekdy_HC51 = TotlPoint_00to12Z_weekdy_file.variables['HC51'][:][:] 
TotlPoint_00to12Z_weekdy_HC52 = TotlPoint_00to12Z_weekdy_file.variables['HC52'][:][:] 
TotlPoint_00to12Z_weekdy_HC53 = TotlPoint_00to12Z_weekdy_file.variables['HC53'][:][:] 
TotlPoint_00to12Z_weekdy_HC54 = TotlPoint_00to12Z_weekdy_file.variables['HC54'][:][:] 
TotlPoint_00to12Z_weekdy_HC55 = TotlPoint_00to12Z_weekdy_file.variables['HC55'][:][:] 
TotlPoint_00to12Z_weekdy_HC56 = TotlPoint_00to12Z_weekdy_file.variables['HC56'][:][:] 
TotlPoint_00to12Z_weekdy_HC57 = TotlPoint_00to12Z_weekdy_file.variables['HC57'][:][:] 
TotlPoint_00to12Z_weekdy_HC58 = TotlPoint_00to12Z_weekdy_file.variables['HC58'][:][:] 
TotlPoint_00to12Z_weekdy_HC59 = TotlPoint_00to12Z_weekdy_file.variables['HC59'][:][:] 
TotlPoint_00to12Z_weekdy_HC60 = TotlPoint_00to12Z_weekdy_file.variables['HC60'][:][:] 
TotlPoint_00to12Z_weekdy_HC61 = TotlPoint_00to12Z_weekdy_file.variables['HC61'][:][:] 
TotlPoint_00to12Z_weekdy_HC62 = TotlPoint_00to12Z_weekdy_file.variables['HC62'][:][:] 
TotlPoint_00to12Z_weekdy_HC63 = TotlPoint_00to12Z_weekdy_file.variables['HC63'][:][:] 
TotlPoint_00to12Z_weekdy_HC64 = TotlPoint_00to12Z_weekdy_file.variables['HC64'][:][:] 
TotlPoint_00to12Z_weekdy_HC65 = TotlPoint_00to12Z_weekdy_file.variables['HC65'][:][:] 
TotlPoint_00to12Z_weekdy_HC66 = TotlPoint_00to12Z_weekdy_file.variables['HC66'][:][:] 
TotlPoint_00to12Z_weekdy_HC67 = TotlPoint_00to12Z_weekdy_file.variables['HC67'][:][:] 
TotlPoint_00to12Z_weekdy_HC68 = TotlPoint_00to12Z_weekdy_file.variables['HC68'][:][:] 
TotlPoint_00to12Z_weekdy_PM01 = TotlPoint_00to12Z_weekdy_file.variables['PM01'][:][:] 
TotlPoint_00to12Z_weekdy_PM02 = TotlPoint_00to12Z_weekdy_file.variables['PM02'][:][:] 
TotlPoint_00to12Z_weekdy_PM03 = TotlPoint_00to12Z_weekdy_file.variables['PM03'][:][:] 
TotlPoint_00to12Z_weekdy_PM04 = TotlPoint_00to12Z_weekdy_file.variables['PM04'][:][:] 
TotlPoint_00to12Z_weekdy_PM05 = TotlPoint_00to12Z_weekdy_file.variables['PM05'][:][:] 
TotlPoint_00to12Z_weekdy_PM06 = TotlPoint_00to12Z_weekdy_file.variables['PM06'][:][:] 
TotlPoint_00to12Z_weekdy_PM07 = TotlPoint_00to12Z_weekdy_file.variables['PM07'][:][:] 
TotlPoint_00to12Z_weekdy_PM08 = TotlPoint_00to12Z_weekdy_file.variables['PM08'][:][:] 
TotlPoint_00to12Z_weekdy_PM09 = TotlPoint_00to12Z_weekdy_file.variables['PM09'][:][:] 
TotlPoint_00to12Z_weekdy_PM10 = TotlPoint_00to12Z_weekdy_file.variables['PM10'][:][:] 
TotlPoint_00to12Z_weekdy_PM11 = TotlPoint_00to12Z_weekdy_file.variables['PM11'][:][:] 
TotlPoint_00to12Z_weekdy_PM12 = TotlPoint_00to12Z_weekdy_file.variables['PM12'][:][:] 
TotlPoint_00to12Z_weekdy_PM13 = TotlPoint_00to12Z_weekdy_file.variables['PM13'][:][:] 
TotlPoint_00to12Z_weekdy_PM14 = TotlPoint_00to12Z_weekdy_file.variables['PM14'][:][:] 
TotlPoint_00to12Z_weekdy_PM15 = TotlPoint_00to12Z_weekdy_file.variables['PM15'][:][:] 
TotlPoint_00to12Z_weekdy_PM16 = TotlPoint_00to12Z_weekdy_file.variables['PM16'][:][:] 
TotlPoint_00to12Z_weekdy_PM17 = TotlPoint_00to12Z_weekdy_file.variables['PM17'][:][:] 
TotlPoint_00to12Z_weekdy_PM18 = TotlPoint_00to12Z_weekdy_file.variables['PM18'][:][:] 
TotlPoint_00to12Z_weekdy_PM19 = TotlPoint_00to12Z_weekdy_file.variables['PM19'][:][:] 
TotlPoint_00to12Z_weekdy_Times = TotlPoint_00to12Z_weekdy_file.variables['Times'][:]

###################################################################################################
#get total ROW
nROW_org, = TotlPoint_00to12Z_weekdy_ITYPE.shape
nROW_extra_OG = len(LON_ng_proc)
nROW_extra = nROW_extra_OG
nROW = nROW_org + nROW_extra
print("nROW_org",nROW_org)
print("nROW_extra",nROW_extra)
print("nROW",nROW)

###################################################################################################
#Organize extra_data
extra_ITYPE_OG = np.array(ERPTYPE_ng_proc)
extra_ITYPE = extra_ITYPE_OG

extra_STKht_OG = np.array(STKHGT_ng_proc)
extra_STKht = extra_STKht_OG

extra_STKdiam_OG = np.array(STKDIAM_ng_proc)
extra_STKdiam = extra_STKdiam_OG

extra_STKtemp_OG = np.array(STKTEMP_ng_proc)
extra_STKtemp = extra_STKtemp_OG

extra_STKve_OG = np.array(STKVEL_ng_proc)
extra_STKve = extra_STKve_OG

extra_STKflw_OG = np.array(STKFLOW_ng_proc)
extra_STKflw = extra_STKflw_OG

extra_FUGht = np.empty(nROW_extra) #FUGht set as empty

extra_XLONG_OG = np.array(LON_ng_proc)
extra_XLONG = extra_XLONG_OG

extra_XLAT_OG = np.array(LAT_ng_proc)
extra_XLAT = extra_XLAT_OG

extra_STATE_OG = STATE_ng_proc

###################################################################################################
#CO2

##################################################################################
extra_CO2_FCPE_ng_proc = HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_weekdy[0:12,:]
extra_CO2_OG = extra_CO2_FCPE_ng_proc

##################################################################################
extra_CO2 =extra_CO2_OG

###################################################################################################
#CH4 from OG can use GHGRP numbers

##################################################################################
extra_CH4_FCPE_ng_proc = HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_weekdy[0:12,:]
extra_CH4_OG = extra_CH4_FCPE_ng_proc

###################################################################################################

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#AQ: using state-level emission ratios to CO2

###################################################################################################
#grab emission ratios from the summary ratio arrays
###################################################################################################

#based on OG process type, species, and state location
extra_X_OG = np.empty([12,nROW_extra_OG,len(species_vector)])
print("extra_X_OG.shape", extra_X_OG.shape)

for pt in range(0,nROW_extra_OG):
    proc_index = 0
    #print("proc_index",proc_index)

    state_abb_cur = extra_STATE_OG[pt]
    
    if state_abb_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_abb_cur)
        #print("state_index",state_index)
    
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * proc_spec_state_emisXdCO2_OG[proc_index,spec_index,state_index]

    else:
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * statistics.mean(proc_spec_state_emisXdCO2_OG[proc_index,spec_index,:])

extra_X_OG_dict = {}
for spec_cur in species_vector:
    #for CH4 use GHGRP numbers
    if spec_cur == 'HC01':
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_CH4_OG
    else:
        spec_index = species_vector.index(spec_cur)
        #print("spec_index",spec_index)
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_X_OG[:,:,spec_index]

###################################################################################################
#stack OG AQ species
extra_X_dict = {}
for spec_cur in species_vector:
    extra_Xi_OG = extra_X_OG_dict["extra_{0}_OG".format(spec_cur)]
    extra_X = extra_Xi_OG
    extra_X_dict["extra_{0}".format(spec_cur)] = extra_X

###################################################################################################
extra_other_spec = np.zeros([12,nROW_extra])

###################################################################################################
#append extra points to original data
TotlPoint_w_extra_00to12Z_weekdy_ITYPE = np.concatenate((TotlPoint_00to12Z_weekdy_ITYPE,extra_ITYPE),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_STKht = np.concatenate((TotlPoint_00to12Z_weekdy_STKht,extra_STKht),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_STKdiam = np.concatenate((TotlPoint_00to12Z_weekdy_STKdiam,extra_STKdiam),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_STKtemp = np.concatenate((TotlPoint_00to12Z_weekdy_STKtemp,extra_STKtemp),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_STKve = np.concatenate((TotlPoint_00to12Z_weekdy_STKve,extra_STKve),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_STKflw = np.concatenate((TotlPoint_00to12Z_weekdy_STKflw,extra_STKflw),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_FUGht = np.concatenate((TotlPoint_00to12Z_weekdy_FUGht,extra_FUGht),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_XLONG = np.concatenate((TotlPoint_00to12Z_weekdy_XLONG,extra_XLONG),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_XLAT = np.concatenate((TotlPoint_00to12Z_weekdy_XLAT,extra_XLAT),axis=0)
TotlPoint_w_extra_00to12Z_weekdy_CO2 = np.concatenate((TotlPoint_00to12Z_weekdy_CO2,extra_CO2),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_CO = np.concatenate((TotlPoint_00to12Z_weekdy_CO,extra_X_dict["extra_CO"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_NH3 = np.concatenate((TotlPoint_00to12Z_weekdy_NH3,extra_X_dict["extra_NH3"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_NOX = np.concatenate((TotlPoint_00to12Z_weekdy_NOX,extra_X_dict["extra_NOX"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM10_PRI = np.concatenate((TotlPoint_00to12Z_weekdy_PM10_PRI,extra_X_dict["extra_PM10-PRI"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM25_PRI = np.concatenate((TotlPoint_00to12Z_weekdy_PM25_PRI,extra_X_dict["extra_PM25-PRI"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_SO2 = np.concatenate((TotlPoint_00to12Z_weekdy_SO2,extra_X_dict["extra_SO2"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_VOC = np.concatenate((TotlPoint_00to12Z_weekdy_VOC,extra_X_dict["extra_VOC"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC01 = np.concatenate((TotlPoint_00to12Z_weekdy_HC01,extra_X_dict["extra_HC01"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC02 = np.concatenate((TotlPoint_00to12Z_weekdy_HC02,extra_X_dict["extra_HC02"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC03 = np.concatenate((TotlPoint_00to12Z_weekdy_HC03,extra_X_dict["extra_HC03"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC04 = np.concatenate((TotlPoint_00to12Z_weekdy_HC04,extra_X_dict["extra_HC04"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC05 = np.concatenate((TotlPoint_00to12Z_weekdy_HC05,extra_X_dict["extra_HC05"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC06 = np.concatenate((TotlPoint_00to12Z_weekdy_HC06,extra_X_dict["extra_HC06"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC07 = np.concatenate((TotlPoint_00to12Z_weekdy_HC07,extra_X_dict["extra_HC07"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC08 = np.concatenate((TotlPoint_00to12Z_weekdy_HC08,extra_X_dict["extra_HC08"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC09 = np.concatenate((TotlPoint_00to12Z_weekdy_HC09,extra_X_dict["extra_HC09"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC10 = np.concatenate((TotlPoint_00to12Z_weekdy_HC10,extra_X_dict["extra_HC10"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC11 = np.concatenate((TotlPoint_00to12Z_weekdy_HC11,extra_X_dict["extra_HC11"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC12 = np.concatenate((TotlPoint_00to12Z_weekdy_HC12,extra_X_dict["extra_HC12"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC13 = np.concatenate((TotlPoint_00to12Z_weekdy_HC13,extra_X_dict["extra_HC13"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC14 = np.concatenate((TotlPoint_00to12Z_weekdy_HC14,extra_X_dict["extra_HC14"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC15 = np.concatenate((TotlPoint_00to12Z_weekdy_HC15,extra_X_dict["extra_HC15"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC16 = np.concatenate((TotlPoint_00to12Z_weekdy_HC16,extra_X_dict["extra_HC16"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC17 = np.concatenate((TotlPoint_00to12Z_weekdy_HC17,extra_X_dict["extra_HC17"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC18 = np.concatenate((TotlPoint_00to12Z_weekdy_HC18,extra_X_dict["extra_HC18"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC19 = np.concatenate((TotlPoint_00to12Z_weekdy_HC19,extra_X_dict["extra_HC19"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC20 = np.concatenate((TotlPoint_00to12Z_weekdy_HC20,extra_X_dict["extra_HC20"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC21 = np.concatenate((TotlPoint_00to12Z_weekdy_HC21,extra_X_dict["extra_HC21"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC22 = np.concatenate((TotlPoint_00to12Z_weekdy_HC22,extra_X_dict["extra_HC22"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC23 = np.concatenate((TotlPoint_00to12Z_weekdy_HC23,extra_X_dict["extra_HC23"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC24 = np.concatenate((TotlPoint_00to12Z_weekdy_HC24,extra_X_dict["extra_HC24"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC25 = np.concatenate((TotlPoint_00to12Z_weekdy_HC25,extra_X_dict["extra_HC25"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC26 = np.concatenate((TotlPoint_00to12Z_weekdy_HC26,extra_X_dict["extra_HC26"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC27 = np.concatenate((TotlPoint_00to12Z_weekdy_HC27,extra_X_dict["extra_HC27"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC28 = np.concatenate((TotlPoint_00to12Z_weekdy_HC28,extra_X_dict["extra_HC28"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC29 = np.concatenate((TotlPoint_00to12Z_weekdy_HC29,extra_X_dict["extra_HC29"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC30 = np.concatenate((TotlPoint_00to12Z_weekdy_HC30,extra_X_dict["extra_HC30"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC31 = np.concatenate((TotlPoint_00to12Z_weekdy_HC31,extra_X_dict["extra_HC31"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC32 = np.concatenate((TotlPoint_00to12Z_weekdy_HC32,extra_X_dict["extra_HC32"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC33 = np.concatenate((TotlPoint_00to12Z_weekdy_HC33,extra_X_dict["extra_HC33"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC34 = np.concatenate((TotlPoint_00to12Z_weekdy_HC34,extra_X_dict["extra_HC34"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC35 = np.concatenate((TotlPoint_00to12Z_weekdy_HC35,extra_X_dict["extra_HC35"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC36 = np.concatenate((TotlPoint_00to12Z_weekdy_HC36,extra_X_dict["extra_HC36"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC37 = np.concatenate((TotlPoint_00to12Z_weekdy_HC37,extra_X_dict["extra_HC37"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC38 = np.concatenate((TotlPoint_00to12Z_weekdy_HC38,extra_X_dict["extra_HC38"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC39 = np.concatenate((TotlPoint_00to12Z_weekdy_HC39,extra_X_dict["extra_HC39"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC40 = np.concatenate((TotlPoint_00to12Z_weekdy_HC40,extra_X_dict["extra_HC40"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC41 = np.concatenate((TotlPoint_00to12Z_weekdy_HC41,extra_X_dict["extra_HC41"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC42 = np.concatenate((TotlPoint_00to12Z_weekdy_HC42,extra_X_dict["extra_HC42"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC43 = np.concatenate((TotlPoint_00to12Z_weekdy_HC43,extra_X_dict["extra_HC43"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC44 = np.concatenate((TotlPoint_00to12Z_weekdy_HC44,extra_X_dict["extra_HC44"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC45 = np.concatenate((TotlPoint_00to12Z_weekdy_HC45,extra_X_dict["extra_HC45"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC46 = np.concatenate((TotlPoint_00to12Z_weekdy_HC46,extra_X_dict["extra_HC46"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC47 = np.concatenate((TotlPoint_00to12Z_weekdy_HC47,extra_X_dict["extra_HC47"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC48 = np.concatenate((TotlPoint_00to12Z_weekdy_HC48,extra_X_dict["extra_HC48"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC49 = np.concatenate((TotlPoint_00to12Z_weekdy_HC49,extra_X_dict["extra_HC49"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC50 = np.concatenate((TotlPoint_00to12Z_weekdy_HC50,extra_X_dict["extra_HC50"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC51 = np.concatenate((TotlPoint_00to12Z_weekdy_HC51,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC52 = np.concatenate((TotlPoint_00to12Z_weekdy_HC52,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC53 = np.concatenate((TotlPoint_00to12Z_weekdy_HC53,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC54 = np.concatenate((TotlPoint_00to12Z_weekdy_HC54,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC55 = np.concatenate((TotlPoint_00to12Z_weekdy_HC55,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC56 = np.concatenate((TotlPoint_00to12Z_weekdy_HC56,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC57 = np.concatenate((TotlPoint_00to12Z_weekdy_HC57,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC58 = np.concatenate((TotlPoint_00to12Z_weekdy_HC58,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC59 = np.concatenate((TotlPoint_00to12Z_weekdy_HC59,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC60 = np.concatenate((TotlPoint_00to12Z_weekdy_HC60,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC61 = np.concatenate((TotlPoint_00to12Z_weekdy_HC61,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC62 = np.concatenate((TotlPoint_00to12Z_weekdy_HC62,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC63 = np.concatenate((TotlPoint_00to12Z_weekdy_HC63,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC64 = np.concatenate((TotlPoint_00to12Z_weekdy_HC64,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC65 = np.concatenate((TotlPoint_00to12Z_weekdy_HC65,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC66 = np.concatenate((TotlPoint_00to12Z_weekdy_HC66,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC67 = np.concatenate((TotlPoint_00to12Z_weekdy_HC67,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_HC68 = np.concatenate((TotlPoint_00to12Z_weekdy_HC68,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM01 = np.concatenate((TotlPoint_00to12Z_weekdy_PM01,extra_X_dict["extra_PM01"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM02 = np.concatenate((TotlPoint_00to12Z_weekdy_PM02,extra_X_dict["extra_PM02"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM03 = np.concatenate((TotlPoint_00to12Z_weekdy_PM03,extra_X_dict["extra_PM03"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM04 = np.concatenate((TotlPoint_00to12Z_weekdy_PM04,extra_X_dict["extra_PM04"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM05 = np.concatenate((TotlPoint_00to12Z_weekdy_PM05,extra_X_dict["extra_PM05"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM06 = np.concatenate((TotlPoint_00to12Z_weekdy_PM06,extra_X_dict["extra_PM06"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM07 = np.concatenate((TotlPoint_00to12Z_weekdy_PM07,extra_X_dict["extra_PM07"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM08 = np.concatenate((TotlPoint_00to12Z_weekdy_PM08,extra_X_dict["extra_PM08"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM09 = np.concatenate((TotlPoint_00to12Z_weekdy_PM09,extra_X_dict["extra_PM09"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM10 = np.concatenate((TotlPoint_00to12Z_weekdy_PM10,extra_X_dict["extra_PM10"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM11 = np.concatenate((TotlPoint_00to12Z_weekdy_PM11,extra_X_dict["extra_PM11"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM12 = np.concatenate((TotlPoint_00to12Z_weekdy_PM12,extra_X_dict["extra_PM12"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM13 = np.concatenate((TotlPoint_00to12Z_weekdy_PM13,extra_X_dict["extra_PM13"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM14 = np.concatenate((TotlPoint_00to12Z_weekdy_PM14,extra_X_dict["extra_PM14"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM15 = np.concatenate((TotlPoint_00to12Z_weekdy_PM15,extra_X_dict["extra_PM15"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM16 = np.concatenate((TotlPoint_00to12Z_weekdy_PM16,extra_X_dict["extra_PM16"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM17 = np.concatenate((TotlPoint_00to12Z_weekdy_PM17,extra_X_dict["extra_PM17"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM18 = np.concatenate((TotlPoint_00to12Z_weekdy_PM18,extra_X_dict["extra_PM18"]),axis=1)
TotlPoint_w_extra_00to12Z_weekdy_PM19 = np.concatenate((TotlPoint_00to12Z_weekdy_PM19,extra_X_dict["extra_PM19"]),axis=1)

###################################################################################################
#write total points with extra points appended
TotlPoint_w_extra_00to12Z_weekdy_fn = append_dir+'/weekdy/PtOG_00to12Z.nc'
TotlPoint_w_extra_00to12Z_weekdy_file = Dataset(TotlPoint_w_extra_00to12Z_weekdy_fn,mode='w',format='NETCDF3_64BIT')

#Creat dimensions
TotlPoint_w_extra_00to12Z_weekdy_file.createDimension("ROW", nROW)
TotlPoint_w_extra_00to12Z_weekdy_file.createDimension("Time", 12)
TotlPoint_w_extra_00to12Z_weekdy_file.sync()

#Create variables
#float ITYPE(ROW) ;
ITYPE = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('ITYPE','f4',('ROW'),fill_value = float(0))
ITYPE[:] = TotlPoint_w_extra_00to12Z_weekdy_ITYPE
varattrs=["FieldType","MemoryOrder","description","units","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['ITYPE'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['ITYPE'], varattr);
        setattr(ITYPE, varattr, varattrVal)

#float STKht(ROW) ;
STKht = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('STKht','f4',('ROW'),fill_value = float(0))
STKht[:] = TotlPoint_w_extra_00to12Z_weekdy_STKht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['STKht'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['STKht'], varattr);
        setattr(STKht, varattr, varattrVal)
        
#float STKdiam(ROW) ;
STKdiam = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('STKdiam','f4',('ROW'),fill_value = float(0))
STKdiam[:] = TotlPoint_w_extra_00to12Z_weekdy_STKdiam
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['STKdiam'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['STKdiam'], varattr);
        setattr(STKdiam, varattr, varattrVal)

#float STKtemp(ROW) ;
STKtemp = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('STKtemp','f4',('ROW'),fill_value = float(0))
STKtemp[:] = TotlPoint_w_extra_00to12Z_weekdy_STKtemp
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['STKtemp'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['STKtemp'], varattr);
        setattr(STKtemp, varattr, varattrVal)
        
#float STKve(ROW) ;
STKve = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('STKve','f4',('ROW'),fill_value = float(0))
STKve[:] = TotlPoint_w_extra_00to12Z_weekdy_STKve
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['STKve'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['STKve'], varattr);
        setattr(STKve, varattr, varattrVal)
        
#float STKflw(ROW) ;
STKflw = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('STKflw','f4',('ROW'),fill_value = float(0))
STKflw[:] = TotlPoint_w_extra_00to12Z_weekdy_STKflw
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['STKflw'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['STKflw'], varattr);
        setattr(STKflw, varattr, varattrVal)
        
#float FUGht(ROW) ;
FUGht = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('FUGht','f4',('ROW'),fill_value = float(0))
FUGht[:] = TotlPoint_w_extra_00to12Z_weekdy_FUGht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['FUGht'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['FUGht'], varattr);
        setattr(FUGht, varattr, varattrVal)
        
#float XLONG(ROW) ;
XLONG = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('XLONG','f4',('ROW'),fill_value = float(0))
XLONG[:] = TotlPoint_w_extra_00to12Z_weekdy_XLONG
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['XLONG'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['XLONG'], varattr);
        setattr(XLONG, varattr, varattrVal)
        
#float XLAT(ROW) ;
XLAT = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('XLAT','f4',('ROW'),fill_value = float(0))
XLAT[:] = TotlPoint_w_extra_00to12Z_weekdy_XLAT
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['XLAT'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['XLAT'], varattr);
        setattr(XLAT, varattr, varattrVal)

#float CO2(Time, ROW) ;
CO2 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('CO2','f4',('Time','ROW'))
CO2[:] = TotlPoint_w_extra_00to12Z_weekdy_CO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['CO2'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['CO2'], varattr);
        setattr(CO2, varattr, varattrVal)
        
#float CO(Time, ROW) ;
CO = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('CO','f4',('Time','ROW'))
CO[:] = TotlPoint_w_extra_00to12Z_weekdy_CO
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['CO'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['CO'], varattr);
        setattr(CO, varattr, varattrVal)
        
#float NH3(Time, ROW) ;
NH3 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('NH3','f4',('Time','ROW'))
NH3[:] = TotlPoint_w_extra_00to12Z_weekdy_NH3
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['NH3'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['NH3'], varattr);
        setattr(NH3, varattr, varattrVal)
        
#float NOX(Time, ROW) ;
NOX = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('NOX','f4',('Time','ROW'))
NOX[:] = TotlPoint_w_extra_00to12Z_weekdy_NOX
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['NOX'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['NOX'], varattr);
        setattr(NOX, varattr, varattrVal)
        
#float PM10-PRI(Time, ROW) ;
PM10_PRI = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM10-PRI','f4',('Time','ROW'))
PM10_PRI[:] = TotlPoint_w_extra_00to12Z_weekdy_PM10_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM10-PRI'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM10-PRI'], varattr);
        setattr(PM10_PRI, varattr, varattrVal)
        
#float PM25-PRI(Time, ROW) ;
PM25_PRI = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM25-PRI','f4',('Time','ROW'))
PM25_PRI[:] = TotlPoint_w_extra_00to12Z_weekdy_PM25_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM25-PRI'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM25-PRI'], varattr);
        setattr(PM25_PRI, varattr, varattrVal)
        
#float SO2(Time, ROW) ;
SO2 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('SO2','f4',('Time','ROW'))
SO2[:] = TotlPoint_w_extra_00to12Z_weekdy_SO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['SO2'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['SO2'], varattr);
        setattr(SO2, varattr, varattrVal)
        
#float VOC(Time, ROW) ;
VOC = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('VOC','f4',('Time','ROW'))
VOC[:] = TotlPoint_w_extra_00to12Z_weekdy_VOC
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['VOC'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['VOC'], varattr);
        setattr(VOC, varattr, varattrVal)
        
#float HC01(Time, ROW) ;
HC01 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC01','f4',('Time','ROW'))
HC01[:] = TotlPoint_w_extra_00to12Z_weekdy_HC01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC01'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC01'], varattr);
        setattr(HC01, varattr, varattrVal)
        
#float HC02(Time, ROW) ;
HC02 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC02','f4',('Time','ROW'))
HC02[:] = TotlPoint_w_extra_00to12Z_weekdy_HC02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC02'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC02'], varattr);
        setattr(HC02, varattr, varattrVal)
        
#float HC03(Time, ROW) ;
HC03 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC03','f4',('Time','ROW'))
HC03[:] = TotlPoint_w_extra_00to12Z_weekdy_HC03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC03'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC03'], varattr);
        setattr(HC03, varattr, varattrVal)
        
#float HC04(Time, ROW) ;
HC04 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC04','f4',('Time','ROW'))
HC04[:] = TotlPoint_w_extra_00to12Z_weekdy_HC04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC04'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC04'], varattr);
        setattr(HC04, varattr, varattrVal)
        
#float HC05(Time, ROW) ;
HC05 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC05','f4',('Time','ROW'))
HC05[:] = TotlPoint_w_extra_00to12Z_weekdy_HC05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC05'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC05'], varattr);
        setattr(HC05, varattr, varattrVal)
        
#float HC06(Time, ROW) ;
HC06 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC06','f4',('Time','ROW'))
HC06[:] = TotlPoint_w_extra_00to12Z_weekdy_HC06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC06'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC06'], varattr);
        setattr(HC06, varattr, varattrVal)
        
#float HC07(Time, ROW) ;
HC07 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC07','f4',('Time','ROW'))
HC07[:] = TotlPoint_w_extra_00to12Z_weekdy_HC07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC07'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC07'], varattr);
        setattr(HC07, varattr, varattrVal)
        
#float HC08(Time, ROW) ;
HC08 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC08','f4',('Time','ROW'))
HC08[:] = TotlPoint_w_extra_00to12Z_weekdy_HC08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC08'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC08'], varattr);
        setattr(HC08, varattr, varattrVal)
        
#float HC09(Time, ROW) ;
HC09 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC09','f4',('Time','ROW'))
HC09[:] = TotlPoint_w_extra_00to12Z_weekdy_HC09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC09'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC09'], varattr);
        setattr(HC09, varattr, varattrVal)
        
#float HC10(Time, ROW) ;
HC10 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC10','f4',('Time','ROW'))
HC10[:] = TotlPoint_w_extra_00to12Z_weekdy_HC10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC10'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC10'], varattr);
        setattr(HC10, varattr, varattrVal)

#float HC11(Time, ROW) ;
HC11 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC11','f4',('Time','ROW'))
HC11[:] = TotlPoint_w_extra_00to12Z_weekdy_HC11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC11'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC11'], varattr);
        setattr(HC11, varattr, varattrVal)
        
#float HC12(Time, ROW) ;
HC12 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC12','f4',('Time','ROW'))
HC12[:] = TotlPoint_w_extra_00to12Z_weekdy_HC12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC12'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC12'], varattr);
        setattr(HC12, varattr, varattrVal)
        
#float HC13(Time, ROW) ;
HC13 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC13','f4',('Time','ROW'))
HC13[:] = TotlPoint_w_extra_00to12Z_weekdy_HC13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC13'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC13'], varattr);
        setattr(HC13, varattr, varattrVal)
        
#float HC14(Time, ROW) ;
HC14 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC14','f4',('Time','ROW'))
HC14[:] = TotlPoint_w_extra_00to12Z_weekdy_HC14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC14'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC14'], varattr);
        setattr(HC14, varattr, varattrVal)
        
#float HC15(Time, ROW) ;
HC15 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC15','f4',('Time','ROW'))
HC15[:] = TotlPoint_w_extra_00to12Z_weekdy_HC15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC15'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC15'], varattr);
        setattr(HC15, varattr, varattrVal)
        
#float HC16(Time, ROW) ;
HC16 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC16','f4',('Time','ROW'))
HC16[:] = TotlPoint_w_extra_00to12Z_weekdy_HC16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC16'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC16'], varattr);
        setattr(HC16, varattr, varattrVal)
        
#float HC17(Time, ROW) ;
HC17 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC17','f4',('Time','ROW'))
HC17[:] = TotlPoint_w_extra_00to12Z_weekdy_HC17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC17'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC17'], varattr);
        setattr(HC17, varattr, varattrVal)
        
#float HC18(Time, ROW) ;
HC18 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC18','f4',('Time','ROW'))
HC18[:] = TotlPoint_w_extra_00to12Z_weekdy_HC18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC18'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC18'], varattr);
        setattr(HC18, varattr, varattrVal)
        
#float HC19(Time, ROW) ;
HC19 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC19','f4',('Time','ROW'))
HC19[:] = TotlPoint_w_extra_00to12Z_weekdy_HC19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC19'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC19'], varattr);
        setattr(HC19, varattr, varattrVal)

#float HC20(Time, ROW) ;
HC20 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC20','f4',('Time','ROW'))
HC20[:] = TotlPoint_w_extra_00to12Z_weekdy_HC20
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC20'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC20'], varattr);
        setattr(HC20, varattr, varattrVal)

#float HC21(Time, ROW) ;
HC21 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC21','f4',('Time','ROW'))
HC21[:] = TotlPoint_w_extra_00to12Z_weekdy_HC21
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC21'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC21'], varattr);
        setattr(HC21, varattr, varattrVal)
        
#float HC22(Time, ROW) ;
HC22 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC22','f4',('Time','ROW'))
HC22[:] = TotlPoint_w_extra_00to12Z_weekdy_HC22
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC22'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC22'], varattr);
        setattr(HC22, varattr, varattrVal)
        
#float HC23(Time, ROW) ;
HC23 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC23','f4',('Time','ROW'))
HC23[:] = TotlPoint_w_extra_00to12Z_weekdy_HC23
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC23'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC23'], varattr);
        setattr(HC23, varattr, varattrVal)
        
#float HC24(Time, ROW) ;
HC24 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC24','f4',('Time','ROW'))
HC24[:] = TotlPoint_w_extra_00to12Z_weekdy_HC24
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC24'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC24'], varattr);
        setattr(HC24, varattr, varattrVal)
        
#float HC25(Time, ROW) ;
HC25 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC25','f4',('Time','ROW'))
HC25[:] = TotlPoint_w_extra_00to12Z_weekdy_HC25
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC25'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC25'], varattr);
        setattr(HC25, varattr, varattrVal)
        
#float HC26(Time, ROW) ;
HC26 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC26','f4',('Time','ROW'))
HC26[:] = TotlPoint_w_extra_00to12Z_weekdy_HC26
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC26'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC26'], varattr);
        setattr(HC26, varattr, varattrVal)
        
#float HC27(Time, ROW) ;
HC27 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC27','f4',('Time','ROW'))
HC27[:] = TotlPoint_w_extra_00to12Z_weekdy_HC27
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC27'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC27'], varattr);
        setattr(HC27, varattr, varattrVal)
        
#float HC28(Time, ROW) ;
HC28 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC28','f4',('Time','ROW'))
HC28[:] = TotlPoint_w_extra_00to12Z_weekdy_HC28
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC28'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC28'], varattr);
        setattr(HC28, varattr, varattrVal)
        
#float HC29(Time, ROW) ;
HC29 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC29','f4',('Time','ROW'))
HC29[:] = TotlPoint_w_extra_00to12Z_weekdy_HC29
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC29'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC29'], varattr);
        setattr(HC29, varattr, varattrVal)

#float HC30(Time, ROW) ;
HC30 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC30','f4',('Time','ROW'))
HC30[:] = TotlPoint_w_extra_00to12Z_weekdy_HC30
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC30'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC30'], varattr);
        setattr(HC30, varattr, varattrVal)

#float HC31(Time, ROW) ;
HC31 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC31','f4',('Time','ROW'))
HC31[:] = TotlPoint_w_extra_00to12Z_weekdy_HC31
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC31'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC31'], varattr);
        setattr(HC31, varattr, varattrVal)
        
#float HC32(Time, ROW) ;
HC32 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC32','f4',('Time','ROW'))
HC32[:] = TotlPoint_w_extra_00to12Z_weekdy_HC32
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC32'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC32'], varattr);
        setattr(HC32, varattr, varattrVal)
        
#float HC33(Time, ROW) ;
HC33 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC33','f4',('Time','ROW'))
HC33[:] = TotlPoint_w_extra_00to12Z_weekdy_HC33
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC33'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC33'], varattr);
        setattr(HC33, varattr, varattrVal)
        
#float HC34(Time, ROW) ;
HC34 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC34','f4',('Time','ROW'))
HC34[:] = TotlPoint_w_extra_00to12Z_weekdy_HC34
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC34'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC34'], varattr);
        setattr(HC34, varattr, varattrVal)
        
#float HC35(Time, ROW) ;
HC35 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC35','f4',('Time','ROW'))
HC35[:] = TotlPoint_w_extra_00to12Z_weekdy_HC35
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC35'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC35'], varattr);
        setattr(HC35, varattr, varattrVal)
        
#float HC36(Time, ROW) ;
HC36 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC36','f4',('Time','ROW'))
HC36[:] = TotlPoint_w_extra_00to12Z_weekdy_HC36
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC36'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC36'], varattr);
        setattr(HC36, varattr, varattrVal)
        
#float HC37(Time, ROW) ;
HC37 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC37','f4',('Time','ROW'))
HC37[:] = TotlPoint_w_extra_00to12Z_weekdy_HC37
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC37'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC37'], varattr);
        setattr(HC37, varattr, varattrVal)
        
#float HC38(Time, ROW) ;
HC38 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC38','f4',('Time','ROW'))
HC38[:] = TotlPoint_w_extra_00to12Z_weekdy_HC38
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC38'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC38'], varattr);
        setattr(HC38, varattr, varattrVal)
        
#float HC39(Time, ROW) ;
HC39 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC39','f4',('Time','ROW'))
HC39[:] = TotlPoint_w_extra_00to12Z_weekdy_HC39
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC39'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC39'], varattr);
        setattr(HC39, varattr, varattrVal)
        
#float HC40(Time, ROW) ;
HC40 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC40','f4',('Time','ROW'))
HC40[:] = TotlPoint_w_extra_00to12Z_weekdy_HC40
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC40'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC40'], varattr);
        setattr(HC40, varattr, varattrVal)

#float HC41(Time, ROW) ;
HC41 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC41','f4',('Time','ROW'))
HC41[:] = TotlPoint_w_extra_00to12Z_weekdy_HC41
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC41'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC41'], varattr);
        setattr(HC41, varattr, varattrVal)
        
#float HC42(Time, ROW) ;
HC42 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC42','f4',('Time','ROW'))
HC42[:] = TotlPoint_w_extra_00to12Z_weekdy_HC42
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC42'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC42'], varattr);
        setattr(HC42, varattr, varattrVal)
        
#float HC43(Time, ROW) ;
HC43 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC43','f4',('Time','ROW'))
HC43[:] = TotlPoint_w_extra_00to12Z_weekdy_HC43
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC43'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC43'], varattr);
        setattr(HC43, varattr, varattrVal)
        
#float HC44(Time, ROW) ;
HC44 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC44','f4',('Time','ROW'))
HC44[:] = TotlPoint_w_extra_00to12Z_weekdy_HC44
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC44'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC44'], varattr);
        setattr(HC44, varattr, varattrVal)
        
#float HC45(Time, ROW) ;
HC45 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC45','f4',('Time','ROW'))
HC45[:] = TotlPoint_w_extra_00to12Z_weekdy_HC45
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC45'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC45'], varattr);
        setattr(HC45, varattr, varattrVal)
        
#float HC46(Time, ROW) ;
HC46 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC46','f4',('Time','ROW'))
HC46[:] = TotlPoint_w_extra_00to12Z_weekdy_HC46
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC46'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC46'], varattr);
        setattr(HC46, varattr, varattrVal)
        
#float HC47(Time, ROW) ;
HC47 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC47','f4',('Time','ROW'))
HC47[:] = TotlPoint_w_extra_00to12Z_weekdy_HC47
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC47'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC47'], varattr);
        setattr(HC47, varattr, varattrVal)
        
#float HC48(Time, ROW) ;
HC48 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC48','f4',('Time','ROW'))
HC48[:] = TotlPoint_w_extra_00to12Z_weekdy_HC48
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC48'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC48'], varattr);
        setattr(HC48, varattr, varattrVal)
        
#float HC49(Time, ROW) ;
HC49 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC49','f4',('Time','ROW'))
HC49[:] = TotlPoint_w_extra_00to12Z_weekdy_HC49
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC49'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC49'], varattr);
        setattr(HC49, varattr, varattrVal)
        
#float HC50(Time, ROW) ;
HC50 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC50','f4',('Time','ROW'))
HC50[:] = TotlPoint_w_extra_00to12Z_weekdy_HC50
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC50'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC50'], varattr);
        setattr(HC50, varattr, varattrVal)

#float HC51(Time, ROW) ;
HC51 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC51','f4',('Time','ROW'))
HC51[:] = TotlPoint_w_extra_00to12Z_weekdy_HC51
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC51'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC51'], varattr);
        setattr(HC51, varattr, varattrVal)
        
#float HC52(Time, ROW) ;
HC52 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC52','f4',('Time','ROW'))
HC52[:] = TotlPoint_w_extra_00to12Z_weekdy_HC52
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC52'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC52'], varattr);
        setattr(HC52, varattr, varattrVal)
        
#float HC53(Time, ROW) ;
HC53 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC53','f4',('Time','ROW'))
HC53[:] = TotlPoint_w_extra_00to12Z_weekdy_HC53
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC53'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC53'], varattr);
        setattr(HC53, varattr, varattrVal)
        
#float HC54(Time, ROW) ;
HC54 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC54','f4',('Time','ROW'))
HC54[:] = TotlPoint_w_extra_00to12Z_weekdy_HC54
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC54'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC54'], varattr);
        setattr(HC54, varattr, varattrVal)
        
#float HC55(Time, ROW) ;
HC55 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC55','f4',('Time','ROW'))
HC55[:] = TotlPoint_w_extra_00to12Z_weekdy_HC55
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC55'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC55'], varattr);
        setattr(HC55, varattr, varattrVal)
        
#float HC56(Time, ROW) ;
HC56 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC56','f4',('Time','ROW'))
HC56[:] = TotlPoint_w_extra_00to12Z_weekdy_HC56
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC56'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC56'], varattr);
        setattr(HC56, varattr, varattrVal)
        
#float HC57(Time, ROW) ;
HC57 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC57','f4',('Time','ROW'))
HC57[:] = TotlPoint_w_extra_00to12Z_weekdy_HC57
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC57'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC57'], varattr);
        setattr(HC57, varattr, varattrVal)
        
#float HC58(Time, ROW) ;
HC58 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC58','f4',('Time','ROW'))
HC58[:] = TotlPoint_w_extra_00to12Z_weekdy_HC58
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC58'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC58'], varattr);
        setattr(HC58, varattr, varattrVal)
        
#float HC59(Time, ROW) ;
HC59 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC59','f4',('Time','ROW'))
HC59[:] = TotlPoint_w_extra_00to12Z_weekdy_HC59
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC59'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC59'], varattr);
        setattr(HC59, varattr, varattrVal)

#float HC60(Time, ROW) ;
HC60 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC60','f4',('Time','ROW'))
HC60[:] = TotlPoint_w_extra_00to12Z_weekdy_HC60
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC60'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC60'], varattr);
        setattr(HC60, varattr, varattrVal)

#float HC61(Time, ROW) ;
HC61 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC61','f4',('Time','ROW'))
HC61[:] = TotlPoint_w_extra_00to12Z_weekdy_HC61
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC61'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC61'], varattr);
        setattr(HC61, varattr, varattrVal)
        
#float HC62(Time, ROW) ;
HC62 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC62','f4',('Time','ROW'))
HC62[:] = TotlPoint_w_extra_00to12Z_weekdy_HC62
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC62'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC62'], varattr);
        setattr(HC62, varattr, varattrVal)
        
#float HC63(Time, ROW) ;
HC63 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC63','f4',('Time','ROW'))
HC63[:] = TotlPoint_w_extra_00to12Z_weekdy_HC63
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC63'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC63'], varattr);
        setattr(HC63, varattr, varattrVal)
        
#float HC64(Time, ROW) ;
HC64 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC64','f4',('Time','ROW'))
HC64[:] = TotlPoint_w_extra_00to12Z_weekdy_HC64
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC64'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC64'], varattr);
        setattr(HC64, varattr, varattrVal)
        
#float HC65(Time, ROW) ;
HC65 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC65','f4',('Time','ROW'))
HC65[:] = TotlPoint_w_extra_00to12Z_weekdy_HC65
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC65'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC65'], varattr);
        setattr(HC65, varattr, varattrVal)
        
#float HC66(Time, ROW) ;
HC66 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC66','f4',('Time','ROW'))
HC66[:] = TotlPoint_w_extra_00to12Z_weekdy_HC66
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC66'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC66'], varattr);
        setattr(HC66, varattr, varattrVal)
        
#float HC67(Time, ROW) ;
HC67 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC67','f4',('Time','ROW'))
HC67[:] = TotlPoint_w_extra_00to12Z_weekdy_HC67
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC67'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC67'], varattr);
        setattr(HC67, varattr, varattrVal)
        
#float HC68(Time, ROW) ;
HC68 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('HC68','f4',('Time','ROW'))
HC68[:] = TotlPoint_w_extra_00to12Z_weekdy_HC68
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['HC68'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['HC68'], varattr);
        setattr(HC68, varattr, varattrVal)

#float PM01(Time, ROW) ;
PM01 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM01','f4',('Time','ROW'))
PM01[:] = TotlPoint_w_extra_00to12Z_weekdy_PM01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM01'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM01'], varattr);
        setattr(PM01, varattr, varattrVal)

#float PM02(Time, ROW) ;
PM02 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM02','f4',('Time','ROW'))
PM02[:] = TotlPoint_w_extra_00to12Z_weekdy_PM02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM02'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM02'], varattr);
        setattr(PM02, varattr, varattrVal)
        
#float PM03(Time, ROW) ;
PM03 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM03','f4',('Time','ROW'))
PM03[:] = TotlPoint_w_extra_00to12Z_weekdy_PM03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM03'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM03'], varattr);
        setattr(PM03, varattr, varattrVal)

#float PM04(Time, ROW) ;
PM04 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM04','f4',('Time','ROW'))
PM04[:] = TotlPoint_w_extra_00to12Z_weekdy_PM04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM04'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM04'], varattr);
        setattr(PM04, varattr, varattrVal)
        
#float PM05(Time, ROW) ;
PM05 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM05','f4',('Time','ROW'))
PM05[:] = TotlPoint_w_extra_00to12Z_weekdy_PM05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM05'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM05'], varattr);
        setattr(PM05, varattr, varattrVal)

#float PM06(Time, ROW) ;
PM06 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM06','f4',('Time','ROW'))
PM06[:] = TotlPoint_w_extra_00to12Z_weekdy_PM06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM06'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM06'], varattr);
        setattr(PM06, varattr, varattrVal)
        
#float PM07(Time, ROW) ;
PM07 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM07','f4',('Time','ROW'))
PM07[:] = TotlPoint_w_extra_00to12Z_weekdy_PM07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM07'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM07'], varattr);
        setattr(PM07, varattr, varattrVal)
        
#float PM08(Time, ROW) ;
PM08 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM08','f4',('Time','ROW'))
PM08[:] = TotlPoint_w_extra_00to12Z_weekdy_PM08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM08'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM08'], varattr);
        setattr(PM08, varattr, varattrVal)
        
#float PM09(Time, ROW) ;
PM09 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM09','f4',('Time','ROW'))
PM09[:] = TotlPoint_w_extra_00to12Z_weekdy_PM09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM09'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM09'], varattr);
        setattr(PM09, varattr, varattrVal)
        
#float PM10(Time, ROW) ;
PM10 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM10','f4',('Time','ROW'))
PM10[:] = TotlPoint_w_extra_00to12Z_weekdy_PM10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM10'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM10'], varattr);
        setattr(PM10, varattr, varattrVal)
        
#float PM11(Time, ROW) ;
PM11 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM11','f4',('Time','ROW'))
PM11[:] = TotlPoint_w_extra_00to12Z_weekdy_PM11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM11'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM11'], varattr);
        setattr(PM11, varattr, varattrVal)

#float PM12(Time, ROW) ;
PM12 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM12','f4',('Time','ROW'))
PM12[:] = TotlPoint_w_extra_00to12Z_weekdy_PM12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM12'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM12'], varattr);
        setattr(PM12, varattr, varattrVal)
        
#float PM13(Time, ROW) ;
PM13 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM13','f4',('Time','ROW'))
PM13[:] = TotlPoint_w_extra_00to12Z_weekdy_PM13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM13'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM13'], varattr);
        setattr(PM13, varattr, varattrVal)

#float PM14(Time, ROW) ;
PM14 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM14','f4',('Time','ROW'))
PM14[:] = TotlPoint_w_extra_00to12Z_weekdy_PM14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM14'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM14'], varattr);
        setattr(PM14, varattr, varattrVal)
        
#float PM15(Time, ROW) ;
PM15 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM15','f4',('Time','ROW'))
PM15[:] = TotlPoint_w_extra_00to12Z_weekdy_PM15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM15'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM15'], varattr);
        setattr(PM15, varattr, varattrVal)

#float PM16(Time, ROW) ;
PM16 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM16','f4',('Time','ROW'))
PM16[:] = TotlPoint_w_extra_00to12Z_weekdy_PM16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM16'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM16'], varattr);
        setattr(PM16, varattr, varattrVal)
        
#float PM17(Time, ROW) ;
PM17 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM17','f4',('Time','ROW'))
PM17[:] = TotlPoint_w_extra_00to12Z_weekdy_PM17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM17'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM17'], varattr);
        setattr(PM17, varattr, varattrVal)
        
#float PM18(Time, ROW) ;
PM18 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM18','f4',('Time','ROW'))
PM18[:] = TotlPoint_w_extra_00to12Z_weekdy_PM18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM18'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM18'], varattr);
        setattr(PM18, varattr, varattrVal)
        
#float PM19(Time, ROW) ;
PM19 = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('PM19','f4',('Time','ROW'))
PM19[:] = TotlPoint_w_extra_00to12Z_weekdy_PM19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_weekdy_file.variables['PM19'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file.variables['PM19'], varattr);
        setattr(PM19, varattr, varattrVal)

#char Times(Time) ;
Times = TotlPoint_w_extra_00to12Z_weekdy_file.createVariable('Times','S1',('Time'))
Times[:] = TotlPoint_00to12Z_weekdy_Times

#copy global attributes from TotlPoint_00to12Z_weekdy_file
for varattr in TotlPoint_00to12Z_weekdy_file.ncattrs():
    if hasattr(TotlPoint_00to12Z_weekdy_file, varattr):
        varattrVal = getattr(TotlPoint_00to12Z_weekdy_file, varattr);
        setattr(TotlPoint_w_extra_00to12Z_weekdy_file, varattr, varattrVal)

TotlPoint_w_extra_00to12Z_weekdy_file.close()


# In[17]:


#append extra points to original point file that is input to wrfchemi assembly program

###################################################################################################
#weekdy, 12to24Z

###################################################################################################
#read original variables
TotlPoint_12to24Z_weekdy_fn = base_dir+'/weekdy/PtOG_12to24Z.nc'
TotlPoint_12to24Z_weekdy_file = Dataset(TotlPoint_12to24Z_weekdy_fn,mode='r',open=True)
TotlPoint_12to24Z_weekdy_ITYPE = TotlPoint_12to24Z_weekdy_file.variables['ITYPE'][:]
TotlPoint_12to24Z_weekdy_STKht = TotlPoint_12to24Z_weekdy_file.variables['STKht'][:]
TotlPoint_12to24Z_weekdy_STKdiam = TotlPoint_12to24Z_weekdy_file.variables['STKdiam'][:]
TotlPoint_12to24Z_weekdy_STKtemp = TotlPoint_12to24Z_weekdy_file.variables['STKtemp'][:]
TotlPoint_12to24Z_weekdy_STKve = TotlPoint_12to24Z_weekdy_file.variables['STKve'][:]
TotlPoint_12to24Z_weekdy_STKflw = TotlPoint_12to24Z_weekdy_file.variables['STKflw'][:]
TotlPoint_12to24Z_weekdy_FUGht = TotlPoint_12to24Z_weekdy_file.variables['FUGht'][:]
TotlPoint_12to24Z_weekdy_XLONG = TotlPoint_12to24Z_weekdy_file.variables['XLONG'][:]
TotlPoint_12to24Z_weekdy_XLAT = TotlPoint_12to24Z_weekdy_file.variables['XLAT'][:]
TotlPoint_12to24Z_weekdy_CO2 = TotlPoint_12to24Z_weekdy_file.variables['CO2'][:][:] 
TotlPoint_12to24Z_weekdy_CO = TotlPoint_12to24Z_weekdy_file.variables['CO'][:][:] 
TotlPoint_12to24Z_weekdy_NH3 = TotlPoint_12to24Z_weekdy_file.variables['NH3'][:][:] 
TotlPoint_12to24Z_weekdy_NOX = TotlPoint_12to24Z_weekdy_file.variables['NOX'][:][:] 
TotlPoint_12to24Z_weekdy_PM10_PRI = TotlPoint_12to24Z_weekdy_file.variables['PM10-PRI'][:][:] 
TotlPoint_12to24Z_weekdy_PM25_PRI = TotlPoint_12to24Z_weekdy_file.variables['PM25-PRI'][:][:] 
TotlPoint_12to24Z_weekdy_SO2 = TotlPoint_12to24Z_weekdy_file.variables['SO2'][:][:] 
TotlPoint_12to24Z_weekdy_VOC = TotlPoint_12to24Z_weekdy_file.variables['VOC'][:][:] 
TotlPoint_12to24Z_weekdy_HC01 = TotlPoint_12to24Z_weekdy_file.variables['HC01'][:][:] 
TotlPoint_12to24Z_weekdy_HC02 = TotlPoint_12to24Z_weekdy_file.variables['HC02'][:][:] 
TotlPoint_12to24Z_weekdy_HC03 = TotlPoint_12to24Z_weekdy_file.variables['HC03'][:][:] 
TotlPoint_12to24Z_weekdy_HC04 = TotlPoint_12to24Z_weekdy_file.variables['HC04'][:][:] 
TotlPoint_12to24Z_weekdy_HC05 = TotlPoint_12to24Z_weekdy_file.variables['HC05'][:][:] 
TotlPoint_12to24Z_weekdy_HC06 = TotlPoint_12to24Z_weekdy_file.variables['HC06'][:][:] 
TotlPoint_12to24Z_weekdy_HC07 = TotlPoint_12to24Z_weekdy_file.variables['HC07'][:][:] 
TotlPoint_12to24Z_weekdy_HC08 = TotlPoint_12to24Z_weekdy_file.variables['HC08'][:][:] 
TotlPoint_12to24Z_weekdy_HC09 = TotlPoint_12to24Z_weekdy_file.variables['HC09'][:][:] 
TotlPoint_12to24Z_weekdy_HC10 = TotlPoint_12to24Z_weekdy_file.variables['HC10'][:][:] 
TotlPoint_12to24Z_weekdy_HC11 = TotlPoint_12to24Z_weekdy_file.variables['HC11'][:][:] 
TotlPoint_12to24Z_weekdy_HC12 = TotlPoint_12to24Z_weekdy_file.variables['HC12'][:][:] 
TotlPoint_12to24Z_weekdy_HC13 = TotlPoint_12to24Z_weekdy_file.variables['HC13'][:][:] 
TotlPoint_12to24Z_weekdy_HC14 = TotlPoint_12to24Z_weekdy_file.variables['HC14'][:][:] 
TotlPoint_12to24Z_weekdy_HC15 = TotlPoint_12to24Z_weekdy_file.variables['HC15'][:][:] 
TotlPoint_12to24Z_weekdy_HC16 = TotlPoint_12to24Z_weekdy_file.variables['HC16'][:][:] 
TotlPoint_12to24Z_weekdy_HC17 = TotlPoint_12to24Z_weekdy_file.variables['HC17'][:][:] 
TotlPoint_12to24Z_weekdy_HC18 = TotlPoint_12to24Z_weekdy_file.variables['HC18'][:][:] 
TotlPoint_12to24Z_weekdy_HC19 = TotlPoint_12to24Z_weekdy_file.variables['HC19'][:][:] 
TotlPoint_12to24Z_weekdy_HC20 = TotlPoint_12to24Z_weekdy_file.variables['HC20'][:][:] 
TotlPoint_12to24Z_weekdy_HC21 = TotlPoint_12to24Z_weekdy_file.variables['HC21'][:][:] 
TotlPoint_12to24Z_weekdy_HC22 = TotlPoint_12to24Z_weekdy_file.variables['HC22'][:][:] 
TotlPoint_12to24Z_weekdy_HC23 = TotlPoint_12to24Z_weekdy_file.variables['HC23'][:][:] 
TotlPoint_12to24Z_weekdy_HC24 = TotlPoint_12to24Z_weekdy_file.variables['HC24'][:][:] 
TotlPoint_12to24Z_weekdy_HC25 = TotlPoint_12to24Z_weekdy_file.variables['HC25'][:][:] 
TotlPoint_12to24Z_weekdy_HC26 = TotlPoint_12to24Z_weekdy_file.variables['HC26'][:][:] 
TotlPoint_12to24Z_weekdy_HC27 = TotlPoint_12to24Z_weekdy_file.variables['HC27'][:][:] 
TotlPoint_12to24Z_weekdy_HC28 = TotlPoint_12to24Z_weekdy_file.variables['HC28'][:][:] 
TotlPoint_12to24Z_weekdy_HC29 = TotlPoint_12to24Z_weekdy_file.variables['HC29'][:][:] 
TotlPoint_12to24Z_weekdy_HC30 = TotlPoint_12to24Z_weekdy_file.variables['HC30'][:][:] 
TotlPoint_12to24Z_weekdy_HC31 = TotlPoint_12to24Z_weekdy_file.variables['HC31'][:][:] 
TotlPoint_12to24Z_weekdy_HC32 = TotlPoint_12to24Z_weekdy_file.variables['HC32'][:][:] 
TotlPoint_12to24Z_weekdy_HC33 = TotlPoint_12to24Z_weekdy_file.variables['HC33'][:][:] 
TotlPoint_12to24Z_weekdy_HC34 = TotlPoint_12to24Z_weekdy_file.variables['HC34'][:][:] 
TotlPoint_12to24Z_weekdy_HC35 = TotlPoint_12to24Z_weekdy_file.variables['HC35'][:][:] 
TotlPoint_12to24Z_weekdy_HC36 = TotlPoint_12to24Z_weekdy_file.variables['HC36'][:][:] 
TotlPoint_12to24Z_weekdy_HC37 = TotlPoint_12to24Z_weekdy_file.variables['HC37'][:][:] 
TotlPoint_12to24Z_weekdy_HC38 = TotlPoint_12to24Z_weekdy_file.variables['HC38'][:][:] 
TotlPoint_12to24Z_weekdy_HC39 = TotlPoint_12to24Z_weekdy_file.variables['HC39'][:][:] 
TotlPoint_12to24Z_weekdy_HC40 = TotlPoint_12to24Z_weekdy_file.variables['HC40'][:][:] 
TotlPoint_12to24Z_weekdy_HC41 = TotlPoint_12to24Z_weekdy_file.variables['HC41'][:][:] 
TotlPoint_12to24Z_weekdy_HC42 = TotlPoint_12to24Z_weekdy_file.variables['HC42'][:][:] 
TotlPoint_12to24Z_weekdy_HC43 = TotlPoint_12to24Z_weekdy_file.variables['HC43'][:][:] 
TotlPoint_12to24Z_weekdy_HC44 = TotlPoint_12to24Z_weekdy_file.variables['HC44'][:][:] 
TotlPoint_12to24Z_weekdy_HC45 = TotlPoint_12to24Z_weekdy_file.variables['HC45'][:][:] 
TotlPoint_12to24Z_weekdy_HC46 = TotlPoint_12to24Z_weekdy_file.variables['HC46'][:][:] 
TotlPoint_12to24Z_weekdy_HC47 = TotlPoint_12to24Z_weekdy_file.variables['HC47'][:][:] 
TotlPoint_12to24Z_weekdy_HC48 = TotlPoint_12to24Z_weekdy_file.variables['HC48'][:][:] 
TotlPoint_12to24Z_weekdy_HC49 = TotlPoint_12to24Z_weekdy_file.variables['HC49'][:][:] 
TotlPoint_12to24Z_weekdy_HC50 = TotlPoint_12to24Z_weekdy_file.variables['HC50'][:][:] 
TotlPoint_12to24Z_weekdy_HC51 = TotlPoint_12to24Z_weekdy_file.variables['HC51'][:][:] 
TotlPoint_12to24Z_weekdy_HC52 = TotlPoint_12to24Z_weekdy_file.variables['HC52'][:][:] 
TotlPoint_12to24Z_weekdy_HC53 = TotlPoint_12to24Z_weekdy_file.variables['HC53'][:][:] 
TotlPoint_12to24Z_weekdy_HC54 = TotlPoint_12to24Z_weekdy_file.variables['HC54'][:][:] 
TotlPoint_12to24Z_weekdy_HC55 = TotlPoint_12to24Z_weekdy_file.variables['HC55'][:][:] 
TotlPoint_12to24Z_weekdy_HC56 = TotlPoint_12to24Z_weekdy_file.variables['HC56'][:][:] 
TotlPoint_12to24Z_weekdy_HC57 = TotlPoint_12to24Z_weekdy_file.variables['HC57'][:][:] 
TotlPoint_12to24Z_weekdy_HC58 = TotlPoint_12to24Z_weekdy_file.variables['HC58'][:][:] 
TotlPoint_12to24Z_weekdy_HC59 = TotlPoint_12to24Z_weekdy_file.variables['HC59'][:][:] 
TotlPoint_12to24Z_weekdy_HC60 = TotlPoint_12to24Z_weekdy_file.variables['HC60'][:][:] 
TotlPoint_12to24Z_weekdy_HC61 = TotlPoint_12to24Z_weekdy_file.variables['HC61'][:][:] 
TotlPoint_12to24Z_weekdy_HC62 = TotlPoint_12to24Z_weekdy_file.variables['HC62'][:][:] 
TotlPoint_12to24Z_weekdy_HC63 = TotlPoint_12to24Z_weekdy_file.variables['HC63'][:][:] 
TotlPoint_12to24Z_weekdy_HC64 = TotlPoint_12to24Z_weekdy_file.variables['HC64'][:][:] 
TotlPoint_12to24Z_weekdy_HC65 = TotlPoint_12to24Z_weekdy_file.variables['HC65'][:][:] 
TotlPoint_12to24Z_weekdy_HC66 = TotlPoint_12to24Z_weekdy_file.variables['HC66'][:][:] 
TotlPoint_12to24Z_weekdy_HC67 = TotlPoint_12to24Z_weekdy_file.variables['HC67'][:][:] 
TotlPoint_12to24Z_weekdy_HC68 = TotlPoint_12to24Z_weekdy_file.variables['HC68'][:][:] 
TotlPoint_12to24Z_weekdy_PM01 = TotlPoint_12to24Z_weekdy_file.variables['PM01'][:][:] 
TotlPoint_12to24Z_weekdy_PM02 = TotlPoint_12to24Z_weekdy_file.variables['PM02'][:][:] 
TotlPoint_12to24Z_weekdy_PM03 = TotlPoint_12to24Z_weekdy_file.variables['PM03'][:][:] 
TotlPoint_12to24Z_weekdy_PM04 = TotlPoint_12to24Z_weekdy_file.variables['PM04'][:][:] 
TotlPoint_12to24Z_weekdy_PM05 = TotlPoint_12to24Z_weekdy_file.variables['PM05'][:][:] 
TotlPoint_12to24Z_weekdy_PM06 = TotlPoint_12to24Z_weekdy_file.variables['PM06'][:][:] 
TotlPoint_12to24Z_weekdy_PM07 = TotlPoint_12to24Z_weekdy_file.variables['PM07'][:][:] 
TotlPoint_12to24Z_weekdy_PM08 = TotlPoint_12to24Z_weekdy_file.variables['PM08'][:][:] 
TotlPoint_12to24Z_weekdy_PM09 = TotlPoint_12to24Z_weekdy_file.variables['PM09'][:][:] 
TotlPoint_12to24Z_weekdy_PM10 = TotlPoint_12to24Z_weekdy_file.variables['PM10'][:][:] 
TotlPoint_12to24Z_weekdy_PM11 = TotlPoint_12to24Z_weekdy_file.variables['PM11'][:][:] 
TotlPoint_12to24Z_weekdy_PM12 = TotlPoint_12to24Z_weekdy_file.variables['PM12'][:][:] 
TotlPoint_12to24Z_weekdy_PM13 = TotlPoint_12to24Z_weekdy_file.variables['PM13'][:][:] 
TotlPoint_12to24Z_weekdy_PM14 = TotlPoint_12to24Z_weekdy_file.variables['PM14'][:][:] 
TotlPoint_12to24Z_weekdy_PM15 = TotlPoint_12to24Z_weekdy_file.variables['PM15'][:][:] 
TotlPoint_12to24Z_weekdy_PM16 = TotlPoint_12to24Z_weekdy_file.variables['PM16'][:][:] 
TotlPoint_12to24Z_weekdy_PM17 = TotlPoint_12to24Z_weekdy_file.variables['PM17'][:][:] 
TotlPoint_12to24Z_weekdy_PM18 = TotlPoint_12to24Z_weekdy_file.variables['PM18'][:][:] 
TotlPoint_12to24Z_weekdy_PM19 = TotlPoint_12to24Z_weekdy_file.variables['PM19'][:][:] 
TotlPoint_12to24Z_weekdy_Times = TotlPoint_12to24Z_weekdy_file.variables['Times'][:]

###################################################################################################
#get total ROW
nROW_org, = TotlPoint_12to24Z_weekdy_ITYPE.shape
nROW_extra_OG = len(LON_ng_proc)
nROW_extra = nROW_extra_OG
nROW = nROW_org + nROW_extra
print("nROW_org",nROW_org)
print("nROW_extra",nROW_extra)
print("nROW",nROW)

###################################################################################################
#Organize extra_data
extra_ITYPE_OG = np.array(ERPTYPE_ng_proc)
extra_ITYPE = extra_ITYPE_OG

extra_STKht_OG = np.array(STKHGT_ng_proc)
extra_STKht = extra_STKht_OG

extra_STKdiam_OG = np.array(STKDIAM_ng_proc)
extra_STKdiam = extra_STKdiam_OG

extra_STKtemp_OG = np.array(STKTEMP_ng_proc)
extra_STKtemp = extra_STKtemp_OG

extra_STKve_OG = np.array(STKVEL_ng_proc)
extra_STKve = extra_STKve_OG

extra_STKflw_OG = np.array(STKFLOW_ng_proc)
extra_STKflw = extra_STKflw_OG

extra_FUGht = np.empty(nROW_extra) #FUGht set as empty

extra_XLONG_OG = np.array(LON_ng_proc)
extra_XLONG = extra_XLONG_OG

extra_XLAT_OG = np.array(LAT_ng_proc)
extra_XLAT = extra_XLAT_OG

extra_STATE_OG = STATE_ng_proc

###################################################################################################
#CO2

##################################################################################
extra_CO2_FCPE_ng_proc = HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_weekdy[12:24,:]
extra_CO2_OG = extra_CO2_FCPE_ng_proc

##################################################################################
extra_CO2 =extra_CO2_OG

###################################################################################################
#CH4 from OG can use GHGRP numbers

##################################################################################
extra_CH4_FCPE_ng_proc = HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_weekdy[12:24,:]
extra_CH4_OG = extra_CH4_FCPE_ng_proc

###################################################################################################

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#AQ: using state-level emission ratios to CO2

###################################################################################################
#grab emission ratios from the summary ratio arrays
###################################################################################################

#based on OG process type, species, and state location
extra_X_OG = np.empty([12,nROW_extra_OG,len(species_vector)])
print("extra_X_OG.shape", extra_X_OG.shape)

for pt in range(0,nROW_extra_OG):
    proc_index = 0
    #print("proc_index",proc_index)

    state_abb_cur = extra_STATE_OG[pt]
    
    if state_abb_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_abb_cur)
        #print("state_index",state_index)
    
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * proc_spec_state_emisXdCO2_OG[proc_index,spec_index,state_index]

    else:
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * statistics.mean(proc_spec_state_emisXdCO2_OG[proc_index,spec_index,:])

extra_X_OG_dict = {}
for spec_cur in species_vector:
    #for CH4 use GHGRP numbers
    if spec_cur == 'HC01':
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_CH4_OG
    else:
        spec_index = species_vector.index(spec_cur)
        #print("spec_index",spec_index)
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_X_OG[:,:,spec_index]

###################################################################################################
#stack OG AQ species
extra_X_dict = {}
for spec_cur in species_vector:
    extra_Xi_OG = extra_X_OG_dict["extra_{0}_OG".format(spec_cur)]
    extra_X = extra_Xi_OG
    extra_X_dict["extra_{0}".format(spec_cur)] = extra_X

###################################################################################################
extra_other_spec = np.zeros([12,nROW_extra])

###################################################################################################
#append extra points to original data
TotlPoint_w_extra_12to24Z_weekdy_ITYPE = np.concatenate((TotlPoint_12to24Z_weekdy_ITYPE,extra_ITYPE),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_STKht = np.concatenate((TotlPoint_12to24Z_weekdy_STKht,extra_STKht),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_STKdiam = np.concatenate((TotlPoint_12to24Z_weekdy_STKdiam,extra_STKdiam),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_STKtemp = np.concatenate((TotlPoint_12to24Z_weekdy_STKtemp,extra_STKtemp),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_STKve = np.concatenate((TotlPoint_12to24Z_weekdy_STKve,extra_STKve),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_STKflw = np.concatenate((TotlPoint_12to24Z_weekdy_STKflw,extra_STKflw),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_FUGht = np.concatenate((TotlPoint_12to24Z_weekdy_FUGht,extra_FUGht),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_XLONG = np.concatenate((TotlPoint_12to24Z_weekdy_XLONG,extra_XLONG),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_XLAT = np.concatenate((TotlPoint_12to24Z_weekdy_XLAT,extra_XLAT),axis=0)
TotlPoint_w_extra_12to24Z_weekdy_CO2 = np.concatenate((TotlPoint_12to24Z_weekdy_CO2,extra_CO2),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_CO = np.concatenate((TotlPoint_12to24Z_weekdy_CO,extra_X_dict["extra_CO"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_NH3 = np.concatenate((TotlPoint_12to24Z_weekdy_NH3,extra_X_dict["extra_NH3"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_NOX = np.concatenate((TotlPoint_12to24Z_weekdy_NOX,extra_X_dict["extra_NOX"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM10_PRI = np.concatenate((TotlPoint_12to24Z_weekdy_PM10_PRI,extra_X_dict["extra_PM10-PRI"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM25_PRI = np.concatenate((TotlPoint_12to24Z_weekdy_PM25_PRI,extra_X_dict["extra_PM25-PRI"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_SO2 = np.concatenate((TotlPoint_12to24Z_weekdy_SO2,extra_X_dict["extra_SO2"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_VOC = np.concatenate((TotlPoint_12to24Z_weekdy_VOC,extra_X_dict["extra_VOC"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC01 = np.concatenate((TotlPoint_12to24Z_weekdy_HC01,extra_X_dict["extra_HC01"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC02 = np.concatenate((TotlPoint_12to24Z_weekdy_HC02,extra_X_dict["extra_HC02"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC03 = np.concatenate((TotlPoint_12to24Z_weekdy_HC03,extra_X_dict["extra_HC03"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC04 = np.concatenate((TotlPoint_12to24Z_weekdy_HC04,extra_X_dict["extra_HC04"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC05 = np.concatenate((TotlPoint_12to24Z_weekdy_HC05,extra_X_dict["extra_HC05"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC06 = np.concatenate((TotlPoint_12to24Z_weekdy_HC06,extra_X_dict["extra_HC06"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC07 = np.concatenate((TotlPoint_12to24Z_weekdy_HC07,extra_X_dict["extra_HC07"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC08 = np.concatenate((TotlPoint_12to24Z_weekdy_HC08,extra_X_dict["extra_HC08"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC09 = np.concatenate((TotlPoint_12to24Z_weekdy_HC09,extra_X_dict["extra_HC09"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC10 = np.concatenate((TotlPoint_12to24Z_weekdy_HC10,extra_X_dict["extra_HC10"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC11 = np.concatenate((TotlPoint_12to24Z_weekdy_HC11,extra_X_dict["extra_HC11"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC12 = np.concatenate((TotlPoint_12to24Z_weekdy_HC12,extra_X_dict["extra_HC12"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC13 = np.concatenate((TotlPoint_12to24Z_weekdy_HC13,extra_X_dict["extra_HC13"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC14 = np.concatenate((TotlPoint_12to24Z_weekdy_HC14,extra_X_dict["extra_HC14"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC15 = np.concatenate((TotlPoint_12to24Z_weekdy_HC15,extra_X_dict["extra_HC15"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC16 = np.concatenate((TotlPoint_12to24Z_weekdy_HC16,extra_X_dict["extra_HC16"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC17 = np.concatenate((TotlPoint_12to24Z_weekdy_HC17,extra_X_dict["extra_HC17"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC18 = np.concatenate((TotlPoint_12to24Z_weekdy_HC18,extra_X_dict["extra_HC18"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC19 = np.concatenate((TotlPoint_12to24Z_weekdy_HC19,extra_X_dict["extra_HC19"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC20 = np.concatenate((TotlPoint_12to24Z_weekdy_HC20,extra_X_dict["extra_HC20"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC21 = np.concatenate((TotlPoint_12to24Z_weekdy_HC21,extra_X_dict["extra_HC21"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC22 = np.concatenate((TotlPoint_12to24Z_weekdy_HC22,extra_X_dict["extra_HC22"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC23 = np.concatenate((TotlPoint_12to24Z_weekdy_HC23,extra_X_dict["extra_HC23"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC24 = np.concatenate((TotlPoint_12to24Z_weekdy_HC24,extra_X_dict["extra_HC24"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC25 = np.concatenate((TotlPoint_12to24Z_weekdy_HC25,extra_X_dict["extra_HC25"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC26 = np.concatenate((TotlPoint_12to24Z_weekdy_HC26,extra_X_dict["extra_HC26"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC27 = np.concatenate((TotlPoint_12to24Z_weekdy_HC27,extra_X_dict["extra_HC27"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC28 = np.concatenate((TotlPoint_12to24Z_weekdy_HC28,extra_X_dict["extra_HC28"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC29 = np.concatenate((TotlPoint_12to24Z_weekdy_HC29,extra_X_dict["extra_HC29"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC30 = np.concatenate((TotlPoint_12to24Z_weekdy_HC30,extra_X_dict["extra_HC30"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC31 = np.concatenate((TotlPoint_12to24Z_weekdy_HC31,extra_X_dict["extra_HC31"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC32 = np.concatenate((TotlPoint_12to24Z_weekdy_HC32,extra_X_dict["extra_HC32"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC33 = np.concatenate((TotlPoint_12to24Z_weekdy_HC33,extra_X_dict["extra_HC33"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC34 = np.concatenate((TotlPoint_12to24Z_weekdy_HC34,extra_X_dict["extra_HC34"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC35 = np.concatenate((TotlPoint_12to24Z_weekdy_HC35,extra_X_dict["extra_HC35"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC36 = np.concatenate((TotlPoint_12to24Z_weekdy_HC36,extra_X_dict["extra_HC36"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC37 = np.concatenate((TotlPoint_12to24Z_weekdy_HC37,extra_X_dict["extra_HC37"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC38 = np.concatenate((TotlPoint_12to24Z_weekdy_HC38,extra_X_dict["extra_HC38"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC39 = np.concatenate((TotlPoint_12to24Z_weekdy_HC39,extra_X_dict["extra_HC39"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC40 = np.concatenate((TotlPoint_12to24Z_weekdy_HC40,extra_X_dict["extra_HC40"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC41 = np.concatenate((TotlPoint_12to24Z_weekdy_HC41,extra_X_dict["extra_HC41"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC42 = np.concatenate((TotlPoint_12to24Z_weekdy_HC42,extra_X_dict["extra_HC42"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC43 = np.concatenate((TotlPoint_12to24Z_weekdy_HC43,extra_X_dict["extra_HC43"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC44 = np.concatenate((TotlPoint_12to24Z_weekdy_HC44,extra_X_dict["extra_HC44"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC45 = np.concatenate((TotlPoint_12to24Z_weekdy_HC45,extra_X_dict["extra_HC45"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC46 = np.concatenate((TotlPoint_12to24Z_weekdy_HC46,extra_X_dict["extra_HC46"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC47 = np.concatenate((TotlPoint_12to24Z_weekdy_HC47,extra_X_dict["extra_HC47"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC48 = np.concatenate((TotlPoint_12to24Z_weekdy_HC48,extra_X_dict["extra_HC48"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC49 = np.concatenate((TotlPoint_12to24Z_weekdy_HC49,extra_X_dict["extra_HC49"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC50 = np.concatenate((TotlPoint_12to24Z_weekdy_HC50,extra_X_dict["extra_HC50"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC51 = np.concatenate((TotlPoint_12to24Z_weekdy_HC51,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC52 = np.concatenate((TotlPoint_12to24Z_weekdy_HC52,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC53 = np.concatenate((TotlPoint_12to24Z_weekdy_HC53,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC54 = np.concatenate((TotlPoint_12to24Z_weekdy_HC54,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC55 = np.concatenate((TotlPoint_12to24Z_weekdy_HC55,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC56 = np.concatenate((TotlPoint_12to24Z_weekdy_HC56,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC57 = np.concatenate((TotlPoint_12to24Z_weekdy_HC57,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC58 = np.concatenate((TotlPoint_12to24Z_weekdy_HC58,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC59 = np.concatenate((TotlPoint_12to24Z_weekdy_HC59,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC60 = np.concatenate((TotlPoint_12to24Z_weekdy_HC60,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC61 = np.concatenate((TotlPoint_12to24Z_weekdy_HC61,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC62 = np.concatenate((TotlPoint_12to24Z_weekdy_HC62,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC63 = np.concatenate((TotlPoint_12to24Z_weekdy_HC63,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC64 = np.concatenate((TotlPoint_12to24Z_weekdy_HC64,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC65 = np.concatenate((TotlPoint_12to24Z_weekdy_HC65,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC66 = np.concatenate((TotlPoint_12to24Z_weekdy_HC66,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC67 = np.concatenate((TotlPoint_12to24Z_weekdy_HC67,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_HC68 = np.concatenate((TotlPoint_12to24Z_weekdy_HC68,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM01 = np.concatenate((TotlPoint_12to24Z_weekdy_PM01,extra_X_dict["extra_PM01"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM02 = np.concatenate((TotlPoint_12to24Z_weekdy_PM02,extra_X_dict["extra_PM02"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM03 = np.concatenate((TotlPoint_12to24Z_weekdy_PM03,extra_X_dict["extra_PM03"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM04 = np.concatenate((TotlPoint_12to24Z_weekdy_PM04,extra_X_dict["extra_PM04"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM05 = np.concatenate((TotlPoint_12to24Z_weekdy_PM05,extra_X_dict["extra_PM05"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM06 = np.concatenate((TotlPoint_12to24Z_weekdy_PM06,extra_X_dict["extra_PM06"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM07 = np.concatenate((TotlPoint_12to24Z_weekdy_PM07,extra_X_dict["extra_PM07"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM08 = np.concatenate((TotlPoint_12to24Z_weekdy_PM08,extra_X_dict["extra_PM08"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM09 = np.concatenate((TotlPoint_12to24Z_weekdy_PM09,extra_X_dict["extra_PM09"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM10 = np.concatenate((TotlPoint_12to24Z_weekdy_PM10,extra_X_dict["extra_PM10"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM11 = np.concatenate((TotlPoint_12to24Z_weekdy_PM11,extra_X_dict["extra_PM11"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM12 = np.concatenate((TotlPoint_12to24Z_weekdy_PM12,extra_X_dict["extra_PM12"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM13 = np.concatenate((TotlPoint_12to24Z_weekdy_PM13,extra_X_dict["extra_PM13"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM14 = np.concatenate((TotlPoint_12to24Z_weekdy_PM14,extra_X_dict["extra_PM14"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM15 = np.concatenate((TotlPoint_12to24Z_weekdy_PM15,extra_X_dict["extra_PM15"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM16 = np.concatenate((TotlPoint_12to24Z_weekdy_PM16,extra_X_dict["extra_PM16"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM17 = np.concatenate((TotlPoint_12to24Z_weekdy_PM17,extra_X_dict["extra_PM17"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM18 = np.concatenate((TotlPoint_12to24Z_weekdy_PM18,extra_X_dict["extra_PM18"]),axis=1)
TotlPoint_w_extra_12to24Z_weekdy_PM19 = np.concatenate((TotlPoint_12to24Z_weekdy_PM19,extra_X_dict["extra_PM19"]),axis=1)

###################################################################################################
#write total points with extra points appended
TotlPoint_w_extra_12to24Z_weekdy_fn = append_dir+'/weekdy/PtOG_12to24Z.nc'
TotlPoint_w_extra_12to24Z_weekdy_file = Dataset(TotlPoint_w_extra_12to24Z_weekdy_fn,mode='w',format='NETCDF3_64BIT')

#Creat dimensions
TotlPoint_w_extra_12to24Z_weekdy_file.createDimension("ROW", nROW)
TotlPoint_w_extra_12to24Z_weekdy_file.createDimension("Time", 12)
TotlPoint_w_extra_12to24Z_weekdy_file.sync()

#Create variables
#float ITYPE(ROW) ;
ITYPE = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('ITYPE','f4',('ROW'),fill_value = float(0))
ITYPE[:] = TotlPoint_w_extra_12to24Z_weekdy_ITYPE
varattrs=["FieldType","MemoryOrder","description","units","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['ITYPE'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['ITYPE'], varattr);
        setattr(ITYPE, varattr, varattrVal)

#float STKht(ROW) ;
STKht = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('STKht','f4',('ROW'),fill_value = float(0))
STKht[:] = TotlPoint_w_extra_12to24Z_weekdy_STKht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['STKht'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['STKht'], varattr);
        setattr(STKht, varattr, varattrVal)
        
#float STKdiam(ROW) ;
STKdiam = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('STKdiam','f4',('ROW'),fill_value = float(0))
STKdiam[:] = TotlPoint_w_extra_12to24Z_weekdy_STKdiam
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['STKdiam'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['STKdiam'], varattr);
        setattr(STKdiam, varattr, varattrVal)

#float STKtemp(ROW) ;
STKtemp = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('STKtemp','f4',('ROW'),fill_value = float(0))
STKtemp[:] = TotlPoint_w_extra_12to24Z_weekdy_STKtemp
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['STKtemp'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['STKtemp'], varattr);
        setattr(STKtemp, varattr, varattrVal)
        
#float STKve(ROW) ;
STKve = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('STKve','f4',('ROW'),fill_value = float(0))
STKve[:] = TotlPoint_w_extra_12to24Z_weekdy_STKve
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['STKve'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['STKve'], varattr);
        setattr(STKve, varattr, varattrVal)
        
#float STKflw(ROW) ;
STKflw = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('STKflw','f4',('ROW'),fill_value = float(0))
STKflw[:] = TotlPoint_w_extra_12to24Z_weekdy_STKflw
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['STKflw'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['STKflw'], varattr);
        setattr(STKflw, varattr, varattrVal)
        
#float FUGht(ROW) ;
FUGht = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('FUGht','f4',('ROW'),fill_value = float(0))
FUGht[:] = TotlPoint_w_extra_12to24Z_weekdy_FUGht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['FUGht'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['FUGht'], varattr);
        setattr(FUGht, varattr, varattrVal)
        
#float XLONG(ROW) ;
XLONG = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('XLONG','f4',('ROW'),fill_value = float(0))
XLONG[:] = TotlPoint_w_extra_12to24Z_weekdy_XLONG
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['XLONG'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['XLONG'], varattr);
        setattr(XLONG, varattr, varattrVal)
        
#float XLAT(ROW) ;
XLAT = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('XLAT','f4',('ROW'),fill_value = float(0))
XLAT[:] = TotlPoint_w_extra_12to24Z_weekdy_XLAT
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['XLAT'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['XLAT'], varattr);
        setattr(XLAT, varattr, varattrVal)

#float CO2(Time, ROW) ;
CO2 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('CO2','f4',('Time','ROW'))
CO2[:] = TotlPoint_w_extra_12to24Z_weekdy_CO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['CO2'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['CO2'], varattr);
        setattr(CO2, varattr, varattrVal)
        
#float CO(Time, ROW) ;
CO = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('CO','f4',('Time','ROW'))
CO[:] = TotlPoint_w_extra_12to24Z_weekdy_CO
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['CO'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['CO'], varattr);
        setattr(CO, varattr, varattrVal)
        
#float NH3(Time, ROW) ;
NH3 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('NH3','f4',('Time','ROW'))
NH3[:] = TotlPoint_w_extra_12to24Z_weekdy_NH3
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['NH3'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['NH3'], varattr);
        setattr(NH3, varattr, varattrVal)
        
#float NOX(Time, ROW) ;
NOX = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('NOX','f4',('Time','ROW'))
NOX[:] = TotlPoint_w_extra_12to24Z_weekdy_NOX
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['NOX'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['NOX'], varattr);
        setattr(NOX, varattr, varattrVal)
        
#float PM10-PRI(Time, ROW) ;
PM10_PRI = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM10-PRI','f4',('Time','ROW'))
PM10_PRI[:] = TotlPoint_w_extra_12to24Z_weekdy_PM10_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM10-PRI'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM10-PRI'], varattr);
        setattr(PM10_PRI, varattr, varattrVal)
        
#float PM25-PRI(Time, ROW) ;
PM25_PRI = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM25-PRI','f4',('Time','ROW'))
PM25_PRI[:] = TotlPoint_w_extra_12to24Z_weekdy_PM25_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM25-PRI'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM25-PRI'], varattr);
        setattr(PM25_PRI, varattr, varattrVal)
        
#float SO2(Time, ROW) ;
SO2 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('SO2','f4',('Time','ROW'))
SO2[:] = TotlPoint_w_extra_12to24Z_weekdy_SO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['SO2'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['SO2'], varattr);
        setattr(SO2, varattr, varattrVal)
        
#float VOC(Time, ROW) ;
VOC = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('VOC','f4',('Time','ROW'))
VOC[:] = TotlPoint_w_extra_12to24Z_weekdy_VOC
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['VOC'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['VOC'], varattr);
        setattr(VOC, varattr, varattrVal)
        
#float HC01(Time, ROW) ;
HC01 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC01','f4',('Time','ROW'))
HC01[:] = TotlPoint_w_extra_12to24Z_weekdy_HC01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC01'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC01'], varattr);
        setattr(HC01, varattr, varattrVal)
        
#float HC02(Time, ROW) ;
HC02 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC02','f4',('Time','ROW'))
HC02[:] = TotlPoint_w_extra_12to24Z_weekdy_HC02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC02'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC02'], varattr);
        setattr(HC02, varattr, varattrVal)
        
#float HC03(Time, ROW) ;
HC03 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC03','f4',('Time','ROW'))
HC03[:] = TotlPoint_w_extra_12to24Z_weekdy_HC03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC03'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC03'], varattr);
        setattr(HC03, varattr, varattrVal)
        
#float HC04(Time, ROW) ;
HC04 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC04','f4',('Time','ROW'))
HC04[:] = TotlPoint_w_extra_12to24Z_weekdy_HC04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC04'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC04'], varattr);
        setattr(HC04, varattr, varattrVal)
        
#float HC05(Time, ROW) ;
HC05 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC05','f4',('Time','ROW'))
HC05[:] = TotlPoint_w_extra_12to24Z_weekdy_HC05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC05'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC05'], varattr);
        setattr(HC05, varattr, varattrVal)
        
#float HC06(Time, ROW) ;
HC06 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC06','f4',('Time','ROW'))
HC06[:] = TotlPoint_w_extra_12to24Z_weekdy_HC06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC06'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC06'], varattr);
        setattr(HC06, varattr, varattrVal)
        
#float HC07(Time, ROW) ;
HC07 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC07','f4',('Time','ROW'))
HC07[:] = TotlPoint_w_extra_12to24Z_weekdy_HC07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC07'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC07'], varattr);
        setattr(HC07, varattr, varattrVal)
        
#float HC08(Time, ROW) ;
HC08 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC08','f4',('Time','ROW'))
HC08[:] = TotlPoint_w_extra_12to24Z_weekdy_HC08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC08'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC08'], varattr);
        setattr(HC08, varattr, varattrVal)
        
#float HC09(Time, ROW) ;
HC09 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC09','f4',('Time','ROW'))
HC09[:] = TotlPoint_w_extra_12to24Z_weekdy_HC09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC09'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC09'], varattr);
        setattr(HC09, varattr, varattrVal)
        
#float HC10(Time, ROW) ;
HC10 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC10','f4',('Time','ROW'))
HC10[:] = TotlPoint_w_extra_12to24Z_weekdy_HC10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC10'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC10'], varattr);
        setattr(HC10, varattr, varattrVal)

#float HC11(Time, ROW) ;
HC11 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC11','f4',('Time','ROW'))
HC11[:] = TotlPoint_w_extra_12to24Z_weekdy_HC11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC11'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC11'], varattr);
        setattr(HC11, varattr, varattrVal)
        
#float HC12(Time, ROW) ;
HC12 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC12','f4',('Time','ROW'))
HC12[:] = TotlPoint_w_extra_12to24Z_weekdy_HC12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC12'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC12'], varattr);
        setattr(HC12, varattr, varattrVal)
        
#float HC13(Time, ROW) ;
HC13 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC13','f4',('Time','ROW'))
HC13[:] = TotlPoint_w_extra_12to24Z_weekdy_HC13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC13'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC13'], varattr);
        setattr(HC13, varattr, varattrVal)
        
#float HC14(Time, ROW) ;
HC14 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC14','f4',('Time','ROW'))
HC14[:] = TotlPoint_w_extra_12to24Z_weekdy_HC14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC14'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC14'], varattr);
        setattr(HC14, varattr, varattrVal)
        
#float HC15(Time, ROW) ;
HC15 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC15','f4',('Time','ROW'))
HC15[:] = TotlPoint_w_extra_12to24Z_weekdy_HC15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC15'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC15'], varattr);
        setattr(HC15, varattr, varattrVal)
        
#float HC16(Time, ROW) ;
HC16 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC16','f4',('Time','ROW'))
HC16[:] = TotlPoint_w_extra_12to24Z_weekdy_HC16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC16'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC16'], varattr);
        setattr(HC16, varattr, varattrVal)
        
#float HC17(Time, ROW) ;
HC17 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC17','f4',('Time','ROW'))
HC17[:] = TotlPoint_w_extra_12to24Z_weekdy_HC17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC17'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC17'], varattr);
        setattr(HC17, varattr, varattrVal)
        
#float HC18(Time, ROW) ;
HC18 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC18','f4',('Time','ROW'))
HC18[:] = TotlPoint_w_extra_12to24Z_weekdy_HC18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC18'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC18'], varattr);
        setattr(HC18, varattr, varattrVal)
        
#float HC19(Time, ROW) ;
HC19 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC19','f4',('Time','ROW'))
HC19[:] = TotlPoint_w_extra_12to24Z_weekdy_HC19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC19'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC19'], varattr);
        setattr(HC19, varattr, varattrVal)

#float HC20(Time, ROW) ;
HC20 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC20','f4',('Time','ROW'))
HC20[:] = TotlPoint_w_extra_12to24Z_weekdy_HC20
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC20'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC20'], varattr);
        setattr(HC20, varattr, varattrVal)

#float HC21(Time, ROW) ;
HC21 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC21','f4',('Time','ROW'))
HC21[:] = TotlPoint_w_extra_12to24Z_weekdy_HC21
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC21'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC21'], varattr);
        setattr(HC21, varattr, varattrVal)
        
#float HC22(Time, ROW) ;
HC22 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC22','f4',('Time','ROW'))
HC22[:] = TotlPoint_w_extra_12to24Z_weekdy_HC22
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC22'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC22'], varattr);
        setattr(HC22, varattr, varattrVal)
        
#float HC23(Time, ROW) ;
HC23 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC23','f4',('Time','ROW'))
HC23[:] = TotlPoint_w_extra_12to24Z_weekdy_HC23
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC23'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC23'], varattr);
        setattr(HC23, varattr, varattrVal)
        
#float HC24(Time, ROW) ;
HC24 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC24','f4',('Time','ROW'))
HC24[:] = TotlPoint_w_extra_12to24Z_weekdy_HC24
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC24'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC24'], varattr);
        setattr(HC24, varattr, varattrVal)
        
#float HC25(Time, ROW) ;
HC25 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC25','f4',('Time','ROW'))
HC25[:] = TotlPoint_w_extra_12to24Z_weekdy_HC25
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC25'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC25'], varattr);
        setattr(HC25, varattr, varattrVal)
        
#float HC26(Time, ROW) ;
HC26 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC26','f4',('Time','ROW'))
HC26[:] = TotlPoint_w_extra_12to24Z_weekdy_HC26
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC26'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC26'], varattr);
        setattr(HC26, varattr, varattrVal)
        
#float HC27(Time, ROW) ;
HC27 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC27','f4',('Time','ROW'))
HC27[:] = TotlPoint_w_extra_12to24Z_weekdy_HC27
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC27'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC27'], varattr);
        setattr(HC27, varattr, varattrVal)
        
#float HC28(Time, ROW) ;
HC28 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC28','f4',('Time','ROW'))
HC28[:] = TotlPoint_w_extra_12to24Z_weekdy_HC28
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC28'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC28'], varattr);
        setattr(HC28, varattr, varattrVal)
        
#float HC29(Time, ROW) ;
HC29 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC29','f4',('Time','ROW'))
HC29[:] = TotlPoint_w_extra_12to24Z_weekdy_HC29
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC29'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC29'], varattr);
        setattr(HC29, varattr, varattrVal)

#float HC30(Time, ROW) ;
HC30 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC30','f4',('Time','ROW'))
HC30[:] = TotlPoint_w_extra_12to24Z_weekdy_HC30
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC30'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC30'], varattr);
        setattr(HC30, varattr, varattrVal)

#float HC31(Time, ROW) ;
HC31 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC31','f4',('Time','ROW'))
HC31[:] = TotlPoint_w_extra_12to24Z_weekdy_HC31
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC31'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC31'], varattr);
        setattr(HC31, varattr, varattrVal)
        
#float HC32(Time, ROW) ;
HC32 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC32','f4',('Time','ROW'))
HC32[:] = TotlPoint_w_extra_12to24Z_weekdy_HC32
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC32'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC32'], varattr);
        setattr(HC32, varattr, varattrVal)
        
#float HC33(Time, ROW) ;
HC33 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC33','f4',('Time','ROW'))
HC33[:] = TotlPoint_w_extra_12to24Z_weekdy_HC33
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC33'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC33'], varattr);
        setattr(HC33, varattr, varattrVal)
        
#float HC34(Time, ROW) ;
HC34 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC34','f4',('Time','ROW'))
HC34[:] = TotlPoint_w_extra_12to24Z_weekdy_HC34
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC34'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC34'], varattr);
        setattr(HC34, varattr, varattrVal)
        
#float HC35(Time, ROW) ;
HC35 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC35','f4',('Time','ROW'))
HC35[:] = TotlPoint_w_extra_12to24Z_weekdy_HC35
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC35'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC35'], varattr);
        setattr(HC35, varattr, varattrVal)
        
#float HC36(Time, ROW) ;
HC36 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC36','f4',('Time','ROW'))
HC36[:] = TotlPoint_w_extra_12to24Z_weekdy_HC36
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC36'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC36'], varattr);
        setattr(HC36, varattr, varattrVal)
        
#float HC37(Time, ROW) ;
HC37 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC37','f4',('Time','ROW'))
HC37[:] = TotlPoint_w_extra_12to24Z_weekdy_HC37
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC37'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC37'], varattr);
        setattr(HC37, varattr, varattrVal)
        
#float HC38(Time, ROW) ;
HC38 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC38','f4',('Time','ROW'))
HC38[:] = TotlPoint_w_extra_12to24Z_weekdy_HC38
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC38'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC38'], varattr);
        setattr(HC38, varattr, varattrVal)
        
#float HC39(Time, ROW) ;
HC39 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC39','f4',('Time','ROW'))
HC39[:] = TotlPoint_w_extra_12to24Z_weekdy_HC39
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC39'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC39'], varattr);
        setattr(HC39, varattr, varattrVal)
        
#float HC40(Time, ROW) ;
HC40 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC40','f4',('Time','ROW'))
HC40[:] = TotlPoint_w_extra_12to24Z_weekdy_HC40
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC40'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC40'], varattr);
        setattr(HC40, varattr, varattrVal)

#float HC41(Time, ROW) ;
HC41 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC41','f4',('Time','ROW'))
HC41[:] = TotlPoint_w_extra_12to24Z_weekdy_HC41
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC41'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC41'], varattr);
        setattr(HC41, varattr, varattrVal)
        
#float HC42(Time, ROW) ;
HC42 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC42','f4',('Time','ROW'))
HC42[:] = TotlPoint_w_extra_12to24Z_weekdy_HC42
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC42'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC42'], varattr);
        setattr(HC42, varattr, varattrVal)
        
#float HC43(Time, ROW) ;
HC43 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC43','f4',('Time','ROW'))
HC43[:] = TotlPoint_w_extra_12to24Z_weekdy_HC43
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC43'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC43'], varattr);
        setattr(HC43, varattr, varattrVal)
        
#float HC44(Time, ROW) ;
HC44 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC44','f4',('Time','ROW'))
HC44[:] = TotlPoint_w_extra_12to24Z_weekdy_HC44
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC44'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC44'], varattr);
        setattr(HC44, varattr, varattrVal)
        
#float HC45(Time, ROW) ;
HC45 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC45','f4',('Time','ROW'))
HC45[:] = TotlPoint_w_extra_12to24Z_weekdy_HC45
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC45'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC45'], varattr);
        setattr(HC45, varattr, varattrVal)
        
#float HC46(Time, ROW) ;
HC46 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC46','f4',('Time','ROW'))
HC46[:] = TotlPoint_w_extra_12to24Z_weekdy_HC46
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC46'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC46'], varattr);
        setattr(HC46, varattr, varattrVal)
        
#float HC47(Time, ROW) ;
HC47 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC47','f4',('Time','ROW'))
HC47[:] = TotlPoint_w_extra_12to24Z_weekdy_HC47
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC47'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC47'], varattr);
        setattr(HC47, varattr, varattrVal)
        
#float HC48(Time, ROW) ;
HC48 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC48','f4',('Time','ROW'))
HC48[:] = TotlPoint_w_extra_12to24Z_weekdy_HC48
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC48'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC48'], varattr);
        setattr(HC48, varattr, varattrVal)
        
#float HC49(Time, ROW) ;
HC49 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC49','f4',('Time','ROW'))
HC49[:] = TotlPoint_w_extra_12to24Z_weekdy_HC49
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC49'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC49'], varattr);
        setattr(HC49, varattr, varattrVal)
        
#float HC50(Time, ROW) ;
HC50 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC50','f4',('Time','ROW'))
HC50[:] = TotlPoint_w_extra_12to24Z_weekdy_HC50
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC50'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC50'], varattr);
        setattr(HC50, varattr, varattrVal)

#float HC51(Time, ROW) ;
HC51 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC51','f4',('Time','ROW'))
HC51[:] = TotlPoint_w_extra_12to24Z_weekdy_HC51
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC51'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC51'], varattr);
        setattr(HC51, varattr, varattrVal)
        
#float HC52(Time, ROW) ;
HC52 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC52','f4',('Time','ROW'))
HC52[:] = TotlPoint_w_extra_12to24Z_weekdy_HC52
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC52'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC52'], varattr);
        setattr(HC52, varattr, varattrVal)
        
#float HC53(Time, ROW) ;
HC53 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC53','f4',('Time','ROW'))
HC53[:] = TotlPoint_w_extra_12to24Z_weekdy_HC53
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC53'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC53'], varattr);
        setattr(HC53, varattr, varattrVal)
        
#float HC54(Time, ROW) ;
HC54 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC54','f4',('Time','ROW'))
HC54[:] = TotlPoint_w_extra_12to24Z_weekdy_HC54
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC54'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC54'], varattr);
        setattr(HC54, varattr, varattrVal)
        
#float HC55(Time, ROW) ;
HC55 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC55','f4',('Time','ROW'))
HC55[:] = TotlPoint_w_extra_12to24Z_weekdy_HC55
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC55'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC55'], varattr);
        setattr(HC55, varattr, varattrVal)
        
#float HC56(Time, ROW) ;
HC56 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC56','f4',('Time','ROW'))
HC56[:] = TotlPoint_w_extra_12to24Z_weekdy_HC56
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC56'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC56'], varattr);
        setattr(HC56, varattr, varattrVal)
        
#float HC57(Time, ROW) ;
HC57 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC57','f4',('Time','ROW'))
HC57[:] = TotlPoint_w_extra_12to24Z_weekdy_HC57
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC57'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC57'], varattr);
        setattr(HC57, varattr, varattrVal)
        
#float HC58(Time, ROW) ;
HC58 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC58','f4',('Time','ROW'))
HC58[:] = TotlPoint_w_extra_12to24Z_weekdy_HC58
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC58'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC58'], varattr);
        setattr(HC58, varattr, varattrVal)
        
#float HC59(Time, ROW) ;
HC59 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC59','f4',('Time','ROW'))
HC59[:] = TotlPoint_w_extra_12to24Z_weekdy_HC59
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC59'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC59'], varattr);
        setattr(HC59, varattr, varattrVal)

#float HC60(Time, ROW) ;
HC60 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC60','f4',('Time','ROW'))
HC60[:] = TotlPoint_w_extra_12to24Z_weekdy_HC60
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC60'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC60'], varattr);
        setattr(HC60, varattr, varattrVal)

#float HC61(Time, ROW) ;
HC61 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC61','f4',('Time','ROW'))
HC61[:] = TotlPoint_w_extra_12to24Z_weekdy_HC61
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC61'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC61'], varattr);
        setattr(HC61, varattr, varattrVal)
        
#float HC62(Time, ROW) ;
HC62 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC62','f4',('Time','ROW'))
HC62[:] = TotlPoint_w_extra_12to24Z_weekdy_HC62
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC62'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC62'], varattr);
        setattr(HC62, varattr, varattrVal)
        
#float HC63(Time, ROW) ;
HC63 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC63','f4',('Time','ROW'))
HC63[:] = TotlPoint_w_extra_12to24Z_weekdy_HC63
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC63'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC63'], varattr);
        setattr(HC63, varattr, varattrVal)
        
#float HC64(Time, ROW) ;
HC64 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC64','f4',('Time','ROW'))
HC64[:] = TotlPoint_w_extra_12to24Z_weekdy_HC64
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC64'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC64'], varattr);
        setattr(HC64, varattr, varattrVal)
        
#float HC65(Time, ROW) ;
HC65 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC65','f4',('Time','ROW'))
HC65[:] = TotlPoint_w_extra_12to24Z_weekdy_HC65
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC65'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC65'], varattr);
        setattr(HC65, varattr, varattrVal)
        
#float HC66(Time, ROW) ;
HC66 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC66','f4',('Time','ROW'))
HC66[:] = TotlPoint_w_extra_12to24Z_weekdy_HC66
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC66'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC66'], varattr);
        setattr(HC66, varattr, varattrVal)
        
#float HC67(Time, ROW) ;
HC67 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC67','f4',('Time','ROW'))
HC67[:] = TotlPoint_w_extra_12to24Z_weekdy_HC67
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC67'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC67'], varattr);
        setattr(HC67, varattr, varattrVal)
        
#float HC68(Time, ROW) ;
HC68 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('HC68','f4',('Time','ROW'))
HC68[:] = TotlPoint_w_extra_12to24Z_weekdy_HC68
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['HC68'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['HC68'], varattr);
        setattr(HC68, varattr, varattrVal)

#float PM01(Time, ROW) ;
PM01 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM01','f4',('Time','ROW'))
PM01[:] = TotlPoint_w_extra_12to24Z_weekdy_PM01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM01'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM01'], varattr);
        setattr(PM01, varattr, varattrVal)

#float PM02(Time, ROW) ;
PM02 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM02','f4',('Time','ROW'))
PM02[:] = TotlPoint_w_extra_12to24Z_weekdy_PM02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM02'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM02'], varattr);
        setattr(PM02, varattr, varattrVal)
        
#float PM03(Time, ROW) ;
PM03 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM03','f4',('Time','ROW'))
PM03[:] = TotlPoint_w_extra_12to24Z_weekdy_PM03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM03'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM03'], varattr);
        setattr(PM03, varattr, varattrVal)

#float PM04(Time, ROW) ;
PM04 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM04','f4',('Time','ROW'))
PM04[:] = TotlPoint_w_extra_12to24Z_weekdy_PM04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM04'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM04'], varattr);
        setattr(PM04, varattr, varattrVal)
        
#float PM05(Time, ROW) ;
PM05 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM05','f4',('Time','ROW'))
PM05[:] = TotlPoint_w_extra_12to24Z_weekdy_PM05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM05'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM05'], varattr);
        setattr(PM05, varattr, varattrVal)

#float PM06(Time, ROW) ;
PM06 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM06','f4',('Time','ROW'))
PM06[:] = TotlPoint_w_extra_12to24Z_weekdy_PM06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM06'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM06'], varattr);
        setattr(PM06, varattr, varattrVal)
        
#float PM07(Time, ROW) ;
PM07 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM07','f4',('Time','ROW'))
PM07[:] = TotlPoint_w_extra_12to24Z_weekdy_PM07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM07'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM07'], varattr);
        setattr(PM07, varattr, varattrVal)
        
#float PM08(Time, ROW) ;
PM08 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM08','f4',('Time','ROW'))
PM08[:] = TotlPoint_w_extra_12to24Z_weekdy_PM08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM08'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM08'], varattr);
        setattr(PM08, varattr, varattrVal)
        
#float PM09(Time, ROW) ;
PM09 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM09','f4',('Time','ROW'))
PM09[:] = TotlPoint_w_extra_12to24Z_weekdy_PM09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM09'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM09'], varattr);
        setattr(PM09, varattr, varattrVal)
        
#float PM10(Time, ROW) ;
PM10 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM10','f4',('Time','ROW'))
PM10[:] = TotlPoint_w_extra_12to24Z_weekdy_PM10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM10'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM10'], varattr);
        setattr(PM10, varattr, varattrVal)
        
#float PM11(Time, ROW) ;
PM11 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM11','f4',('Time','ROW'))
PM11[:] = TotlPoint_w_extra_12to24Z_weekdy_PM11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM11'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM11'], varattr);
        setattr(PM11, varattr, varattrVal)

#float PM12(Time, ROW) ;
PM12 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM12','f4',('Time','ROW'))
PM12[:] = TotlPoint_w_extra_12to24Z_weekdy_PM12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM12'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM12'], varattr);
        setattr(PM12, varattr, varattrVal)
        
#float PM13(Time, ROW) ;
PM13 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM13','f4',('Time','ROW'))
PM13[:] = TotlPoint_w_extra_12to24Z_weekdy_PM13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM13'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM13'], varattr);
        setattr(PM13, varattr, varattrVal)

#float PM14(Time, ROW) ;
PM14 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM14','f4',('Time','ROW'))
PM14[:] = TotlPoint_w_extra_12to24Z_weekdy_PM14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM14'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM14'], varattr);
        setattr(PM14, varattr, varattrVal)
        
#float PM15(Time, ROW) ;
PM15 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM15','f4',('Time','ROW'))
PM15[:] = TotlPoint_w_extra_12to24Z_weekdy_PM15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM15'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM15'], varattr);
        setattr(PM15, varattr, varattrVal)

#float PM16(Time, ROW) ;
PM16 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM16','f4',('Time','ROW'))
PM16[:] = TotlPoint_w_extra_12to24Z_weekdy_PM16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM16'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM16'], varattr);
        setattr(PM16, varattr, varattrVal)
        
#float PM17(Time, ROW) ;
PM17 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM17','f4',('Time','ROW'))
PM17[:] = TotlPoint_w_extra_12to24Z_weekdy_PM17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM17'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM17'], varattr);
        setattr(PM17, varattr, varattrVal)
        
#float PM18(Time, ROW) ;
PM18 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM18','f4',('Time','ROW'))
PM18[:] = TotlPoint_w_extra_12to24Z_weekdy_PM18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM18'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM18'], varattr);
        setattr(PM18, varattr, varattrVal)
        
#float PM19(Time, ROW) ;
PM19 = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('PM19','f4',('Time','ROW'))
PM19[:] = TotlPoint_w_extra_12to24Z_weekdy_PM19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_weekdy_file.variables['PM19'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file.variables['PM19'], varattr);
        setattr(PM19, varattr, varattrVal)

#char Times(Time) ;
Times = TotlPoint_w_extra_12to24Z_weekdy_file.createVariable('Times','S1',('Time'))
Times[:] = TotlPoint_12to24Z_weekdy_Times

#copy global attributes from TotlPoint_12to24Z_weekdy_file
for varattr in TotlPoint_12to24Z_weekdy_file.ncattrs():
    if hasattr(TotlPoint_12to24Z_weekdy_file, varattr):
        varattrVal = getattr(TotlPoint_12to24Z_weekdy_file, varattr);
        setattr(TotlPoint_w_extra_12to24Z_weekdy_file, varattr, varattrVal)

TotlPoint_w_extra_12to24Z_weekdy_file.close()


# In[18]:


#append extra points to original point file that is input to wrfchemi assembly program

###################################################################################################
#satdy, 00to12Z

###################################################################################################
#read original variables
TotlPoint_00to12Z_satdy_fn = base_dir+'/satdy/PtOG_00to12Z.nc'
TotlPoint_00to12Z_satdy_file = Dataset(TotlPoint_00to12Z_satdy_fn,mode='r',open=True)
TotlPoint_00to12Z_satdy_ITYPE = TotlPoint_00to12Z_satdy_file.variables['ITYPE'][:]
TotlPoint_00to12Z_satdy_STKht = TotlPoint_00to12Z_satdy_file.variables['STKht'][:]
TotlPoint_00to12Z_satdy_STKdiam = TotlPoint_00to12Z_satdy_file.variables['STKdiam'][:]
TotlPoint_00to12Z_satdy_STKtemp = TotlPoint_00to12Z_satdy_file.variables['STKtemp'][:]
TotlPoint_00to12Z_satdy_STKve = TotlPoint_00to12Z_satdy_file.variables['STKve'][:]
TotlPoint_00to12Z_satdy_STKflw = TotlPoint_00to12Z_satdy_file.variables['STKflw'][:]
TotlPoint_00to12Z_satdy_FUGht = TotlPoint_00to12Z_satdy_file.variables['FUGht'][:]
TotlPoint_00to12Z_satdy_XLONG = TotlPoint_00to12Z_satdy_file.variables['XLONG'][:]
TotlPoint_00to12Z_satdy_XLAT = TotlPoint_00to12Z_satdy_file.variables['XLAT'][:]
TotlPoint_00to12Z_satdy_CO2 = TotlPoint_00to12Z_satdy_file.variables['CO2'][:][:] 
TotlPoint_00to12Z_satdy_CO = TotlPoint_00to12Z_satdy_file.variables['CO'][:][:] 
TotlPoint_00to12Z_satdy_NH3 = TotlPoint_00to12Z_satdy_file.variables['NH3'][:][:] 
TotlPoint_00to12Z_satdy_NOX = TotlPoint_00to12Z_satdy_file.variables['NOX'][:][:] 
TotlPoint_00to12Z_satdy_PM10_PRI = TotlPoint_00to12Z_satdy_file.variables['PM10-PRI'][:][:] 
TotlPoint_00to12Z_satdy_PM25_PRI = TotlPoint_00to12Z_satdy_file.variables['PM25-PRI'][:][:] 
TotlPoint_00to12Z_satdy_SO2 = TotlPoint_00to12Z_satdy_file.variables['SO2'][:][:] 
TotlPoint_00to12Z_satdy_VOC = TotlPoint_00to12Z_satdy_file.variables['VOC'][:][:] 
TotlPoint_00to12Z_satdy_HC01 = TotlPoint_00to12Z_satdy_file.variables['HC01'][:][:] 
TotlPoint_00to12Z_satdy_HC02 = TotlPoint_00to12Z_satdy_file.variables['HC02'][:][:] 
TotlPoint_00to12Z_satdy_HC03 = TotlPoint_00to12Z_satdy_file.variables['HC03'][:][:] 
TotlPoint_00to12Z_satdy_HC04 = TotlPoint_00to12Z_satdy_file.variables['HC04'][:][:] 
TotlPoint_00to12Z_satdy_HC05 = TotlPoint_00to12Z_satdy_file.variables['HC05'][:][:] 
TotlPoint_00to12Z_satdy_HC06 = TotlPoint_00to12Z_satdy_file.variables['HC06'][:][:] 
TotlPoint_00to12Z_satdy_HC07 = TotlPoint_00to12Z_satdy_file.variables['HC07'][:][:] 
TotlPoint_00to12Z_satdy_HC08 = TotlPoint_00to12Z_satdy_file.variables['HC08'][:][:] 
TotlPoint_00to12Z_satdy_HC09 = TotlPoint_00to12Z_satdy_file.variables['HC09'][:][:] 
TotlPoint_00to12Z_satdy_HC10 = TotlPoint_00to12Z_satdy_file.variables['HC10'][:][:] 
TotlPoint_00to12Z_satdy_HC11 = TotlPoint_00to12Z_satdy_file.variables['HC11'][:][:] 
TotlPoint_00to12Z_satdy_HC12 = TotlPoint_00to12Z_satdy_file.variables['HC12'][:][:] 
TotlPoint_00to12Z_satdy_HC13 = TotlPoint_00to12Z_satdy_file.variables['HC13'][:][:] 
TotlPoint_00to12Z_satdy_HC14 = TotlPoint_00to12Z_satdy_file.variables['HC14'][:][:] 
TotlPoint_00to12Z_satdy_HC15 = TotlPoint_00to12Z_satdy_file.variables['HC15'][:][:] 
TotlPoint_00to12Z_satdy_HC16 = TotlPoint_00to12Z_satdy_file.variables['HC16'][:][:] 
TotlPoint_00to12Z_satdy_HC17 = TotlPoint_00to12Z_satdy_file.variables['HC17'][:][:] 
TotlPoint_00to12Z_satdy_HC18 = TotlPoint_00to12Z_satdy_file.variables['HC18'][:][:] 
TotlPoint_00to12Z_satdy_HC19 = TotlPoint_00to12Z_satdy_file.variables['HC19'][:][:] 
TotlPoint_00to12Z_satdy_HC20 = TotlPoint_00to12Z_satdy_file.variables['HC20'][:][:] 
TotlPoint_00to12Z_satdy_HC21 = TotlPoint_00to12Z_satdy_file.variables['HC21'][:][:] 
TotlPoint_00to12Z_satdy_HC22 = TotlPoint_00to12Z_satdy_file.variables['HC22'][:][:] 
TotlPoint_00to12Z_satdy_HC23 = TotlPoint_00to12Z_satdy_file.variables['HC23'][:][:] 
TotlPoint_00to12Z_satdy_HC24 = TotlPoint_00to12Z_satdy_file.variables['HC24'][:][:] 
TotlPoint_00to12Z_satdy_HC25 = TotlPoint_00to12Z_satdy_file.variables['HC25'][:][:] 
TotlPoint_00to12Z_satdy_HC26 = TotlPoint_00to12Z_satdy_file.variables['HC26'][:][:] 
TotlPoint_00to12Z_satdy_HC27 = TotlPoint_00to12Z_satdy_file.variables['HC27'][:][:] 
TotlPoint_00to12Z_satdy_HC28 = TotlPoint_00to12Z_satdy_file.variables['HC28'][:][:] 
TotlPoint_00to12Z_satdy_HC29 = TotlPoint_00to12Z_satdy_file.variables['HC29'][:][:] 
TotlPoint_00to12Z_satdy_HC30 = TotlPoint_00to12Z_satdy_file.variables['HC30'][:][:] 
TotlPoint_00to12Z_satdy_HC31 = TotlPoint_00to12Z_satdy_file.variables['HC31'][:][:] 
TotlPoint_00to12Z_satdy_HC32 = TotlPoint_00to12Z_satdy_file.variables['HC32'][:][:] 
TotlPoint_00to12Z_satdy_HC33 = TotlPoint_00to12Z_satdy_file.variables['HC33'][:][:] 
TotlPoint_00to12Z_satdy_HC34 = TotlPoint_00to12Z_satdy_file.variables['HC34'][:][:] 
TotlPoint_00to12Z_satdy_HC35 = TotlPoint_00to12Z_satdy_file.variables['HC35'][:][:] 
TotlPoint_00to12Z_satdy_HC36 = TotlPoint_00to12Z_satdy_file.variables['HC36'][:][:] 
TotlPoint_00to12Z_satdy_HC37 = TotlPoint_00to12Z_satdy_file.variables['HC37'][:][:] 
TotlPoint_00to12Z_satdy_HC38 = TotlPoint_00to12Z_satdy_file.variables['HC38'][:][:] 
TotlPoint_00to12Z_satdy_HC39 = TotlPoint_00to12Z_satdy_file.variables['HC39'][:][:] 
TotlPoint_00to12Z_satdy_HC40 = TotlPoint_00to12Z_satdy_file.variables['HC40'][:][:] 
TotlPoint_00to12Z_satdy_HC41 = TotlPoint_00to12Z_satdy_file.variables['HC41'][:][:] 
TotlPoint_00to12Z_satdy_HC42 = TotlPoint_00to12Z_satdy_file.variables['HC42'][:][:] 
TotlPoint_00to12Z_satdy_HC43 = TotlPoint_00to12Z_satdy_file.variables['HC43'][:][:] 
TotlPoint_00to12Z_satdy_HC44 = TotlPoint_00to12Z_satdy_file.variables['HC44'][:][:] 
TotlPoint_00to12Z_satdy_HC45 = TotlPoint_00to12Z_satdy_file.variables['HC45'][:][:] 
TotlPoint_00to12Z_satdy_HC46 = TotlPoint_00to12Z_satdy_file.variables['HC46'][:][:] 
TotlPoint_00to12Z_satdy_HC47 = TotlPoint_00to12Z_satdy_file.variables['HC47'][:][:] 
TotlPoint_00to12Z_satdy_HC48 = TotlPoint_00to12Z_satdy_file.variables['HC48'][:][:] 
TotlPoint_00to12Z_satdy_HC49 = TotlPoint_00to12Z_satdy_file.variables['HC49'][:][:] 
TotlPoint_00to12Z_satdy_HC50 = TotlPoint_00to12Z_satdy_file.variables['HC50'][:][:] 
TotlPoint_00to12Z_satdy_HC51 = TotlPoint_00to12Z_satdy_file.variables['HC51'][:][:] 
TotlPoint_00to12Z_satdy_HC52 = TotlPoint_00to12Z_satdy_file.variables['HC52'][:][:] 
TotlPoint_00to12Z_satdy_HC53 = TotlPoint_00to12Z_satdy_file.variables['HC53'][:][:] 
TotlPoint_00to12Z_satdy_HC54 = TotlPoint_00to12Z_satdy_file.variables['HC54'][:][:] 
TotlPoint_00to12Z_satdy_HC55 = TotlPoint_00to12Z_satdy_file.variables['HC55'][:][:] 
TotlPoint_00to12Z_satdy_HC56 = TotlPoint_00to12Z_satdy_file.variables['HC56'][:][:] 
TotlPoint_00to12Z_satdy_HC57 = TotlPoint_00to12Z_satdy_file.variables['HC57'][:][:] 
TotlPoint_00to12Z_satdy_HC58 = TotlPoint_00to12Z_satdy_file.variables['HC58'][:][:] 
TotlPoint_00to12Z_satdy_HC59 = TotlPoint_00to12Z_satdy_file.variables['HC59'][:][:] 
TotlPoint_00to12Z_satdy_HC60 = TotlPoint_00to12Z_satdy_file.variables['HC60'][:][:] 
TotlPoint_00to12Z_satdy_HC61 = TotlPoint_00to12Z_satdy_file.variables['HC61'][:][:] 
TotlPoint_00to12Z_satdy_HC62 = TotlPoint_00to12Z_satdy_file.variables['HC62'][:][:] 
TotlPoint_00to12Z_satdy_HC63 = TotlPoint_00to12Z_satdy_file.variables['HC63'][:][:] 
TotlPoint_00to12Z_satdy_HC64 = TotlPoint_00to12Z_satdy_file.variables['HC64'][:][:] 
TotlPoint_00to12Z_satdy_HC65 = TotlPoint_00to12Z_satdy_file.variables['HC65'][:][:] 
TotlPoint_00to12Z_satdy_HC66 = TotlPoint_00to12Z_satdy_file.variables['HC66'][:][:] 
TotlPoint_00to12Z_satdy_HC67 = TotlPoint_00to12Z_satdy_file.variables['HC67'][:][:] 
TotlPoint_00to12Z_satdy_HC68 = TotlPoint_00to12Z_satdy_file.variables['HC68'][:][:] 
TotlPoint_00to12Z_satdy_PM01 = TotlPoint_00to12Z_satdy_file.variables['PM01'][:][:] 
TotlPoint_00to12Z_satdy_PM02 = TotlPoint_00to12Z_satdy_file.variables['PM02'][:][:] 
TotlPoint_00to12Z_satdy_PM03 = TotlPoint_00to12Z_satdy_file.variables['PM03'][:][:] 
TotlPoint_00to12Z_satdy_PM04 = TotlPoint_00to12Z_satdy_file.variables['PM04'][:][:] 
TotlPoint_00to12Z_satdy_PM05 = TotlPoint_00to12Z_satdy_file.variables['PM05'][:][:] 
TotlPoint_00to12Z_satdy_PM06 = TotlPoint_00to12Z_satdy_file.variables['PM06'][:][:] 
TotlPoint_00to12Z_satdy_PM07 = TotlPoint_00to12Z_satdy_file.variables['PM07'][:][:] 
TotlPoint_00to12Z_satdy_PM08 = TotlPoint_00to12Z_satdy_file.variables['PM08'][:][:] 
TotlPoint_00to12Z_satdy_PM09 = TotlPoint_00to12Z_satdy_file.variables['PM09'][:][:] 
TotlPoint_00to12Z_satdy_PM10 = TotlPoint_00to12Z_satdy_file.variables['PM10'][:][:] 
TotlPoint_00to12Z_satdy_PM11 = TotlPoint_00to12Z_satdy_file.variables['PM11'][:][:] 
TotlPoint_00to12Z_satdy_PM12 = TotlPoint_00to12Z_satdy_file.variables['PM12'][:][:] 
TotlPoint_00to12Z_satdy_PM13 = TotlPoint_00to12Z_satdy_file.variables['PM13'][:][:] 
TotlPoint_00to12Z_satdy_PM14 = TotlPoint_00to12Z_satdy_file.variables['PM14'][:][:] 
TotlPoint_00to12Z_satdy_PM15 = TotlPoint_00to12Z_satdy_file.variables['PM15'][:][:] 
TotlPoint_00to12Z_satdy_PM16 = TotlPoint_00to12Z_satdy_file.variables['PM16'][:][:] 
TotlPoint_00to12Z_satdy_PM17 = TotlPoint_00to12Z_satdy_file.variables['PM17'][:][:] 
TotlPoint_00to12Z_satdy_PM18 = TotlPoint_00to12Z_satdy_file.variables['PM18'][:][:] 
TotlPoint_00to12Z_satdy_PM19 = TotlPoint_00to12Z_satdy_file.variables['PM19'][:][:] 
TotlPoint_00to12Z_satdy_Times = TotlPoint_00to12Z_satdy_file.variables['Times'][:]

###################################################################################################
#get total ROW
nROW_org, = TotlPoint_00to12Z_satdy_ITYPE.shape
nROW_extra_OG = len(LON_ng_proc)
nROW_extra = nROW_extra_OG
nROW = nROW_org + nROW_extra
print("nROW_org",nROW_org)
print("nROW_extra",nROW_extra)
print("nROW",nROW)

###################################################################################################
#Organize extra_data
extra_ITYPE_OG = np.array(ERPTYPE_ng_proc)
extra_ITYPE = extra_ITYPE_OG

extra_STKht_OG = np.array(STKHGT_ng_proc)
extra_STKht = extra_STKht_OG

extra_STKdiam_OG = np.array(STKDIAM_ng_proc)
extra_STKdiam = extra_STKdiam_OG

extra_STKtemp_OG = np.array(STKTEMP_ng_proc)
extra_STKtemp = extra_STKtemp_OG

extra_STKve_OG = np.array(STKVEL_ng_proc)
extra_STKve = extra_STKve_OG

extra_STKflw_OG = np.array(STKFLOW_ng_proc)
extra_STKflw = extra_STKflw_OG

extra_FUGht = np.empty(nROW_extra) #FUGht set as empty

extra_XLONG_OG = np.array(LON_ng_proc)
extra_XLONG = extra_XLONG_OG

extra_XLAT_OG = np.array(LAT_ng_proc)
extra_XLAT = extra_XLAT_OG

extra_STATE_OG = STATE_ng_proc

###################################################################################################
#CO2

##################################################################################
extra_CO2_FCPE_ng_proc = HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_satdy[0:12,:]
extra_CO2_OG = extra_CO2_FCPE_ng_proc

##################################################################################
extra_CO2 =extra_CO2_OG

###################################################################################################
#CH4 from OG can use GHGRP numbers

##################################################################################
extra_CH4_FCPE_ng_proc = HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_satdy[0:12,:]
extra_CH4_OG = extra_CH4_FCPE_ng_proc

###################################################################################################

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#AQ: using state-level emission ratios to CO2

###################################################################################################
#grab emission ratios from the summary ratio arrays
###################################################################################################

#based on OG process type, species, and state location
extra_X_OG = np.empty([12,nROW_extra_OG,len(species_vector)])
print("extra_X_OG.shape", extra_X_OG.shape)

for pt in range(0,nROW_extra_OG):
    proc_index = 0
    #print("proc_index",proc_index)

    state_abb_cur = extra_STATE_OG[pt]
    
    if state_abb_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_abb_cur)
        #print("state_index",state_index)
    
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * proc_spec_state_emisXdCO2_OG[proc_index,spec_index,state_index]

    else:
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * statistics.mean(proc_spec_state_emisXdCO2_OG[proc_index,spec_index,:])

extra_X_OG_dict = {}
for spec_cur in species_vector:
    #for CH4 use GHGRP numbers
    if spec_cur == 'HC01':
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_CH4_OG
    else:
        spec_index = species_vector.index(spec_cur)
        #print("spec_index",spec_index)
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_X_OG[:,:,spec_index]

###################################################################################################
#stack OG AQ species
extra_X_dict = {}
for spec_cur in species_vector:
    extra_Xi_OG = extra_X_OG_dict["extra_{0}_OG".format(spec_cur)]
    extra_X = extra_Xi_OG
    extra_X_dict["extra_{0}".format(spec_cur)] = extra_X

###################################################################################################
extra_other_spec = np.zeros([12,nROW_extra])

###################################################################################################
#append extra points to original data
TotlPoint_w_extra_00to12Z_satdy_ITYPE = np.concatenate((TotlPoint_00to12Z_satdy_ITYPE,extra_ITYPE),axis=0)
TotlPoint_w_extra_00to12Z_satdy_STKht = np.concatenate((TotlPoint_00to12Z_satdy_STKht,extra_STKht),axis=0)
TotlPoint_w_extra_00to12Z_satdy_STKdiam = np.concatenate((TotlPoint_00to12Z_satdy_STKdiam,extra_STKdiam),axis=0)
TotlPoint_w_extra_00to12Z_satdy_STKtemp = np.concatenate((TotlPoint_00to12Z_satdy_STKtemp,extra_STKtemp),axis=0)
TotlPoint_w_extra_00to12Z_satdy_STKve = np.concatenate((TotlPoint_00to12Z_satdy_STKve,extra_STKve),axis=0)
TotlPoint_w_extra_00to12Z_satdy_STKflw = np.concatenate((TotlPoint_00to12Z_satdy_STKflw,extra_STKflw),axis=0)
TotlPoint_w_extra_00to12Z_satdy_FUGht = np.concatenate((TotlPoint_00to12Z_satdy_FUGht,extra_FUGht),axis=0)
TotlPoint_w_extra_00to12Z_satdy_XLONG = np.concatenate((TotlPoint_00to12Z_satdy_XLONG,extra_XLONG),axis=0)
TotlPoint_w_extra_00to12Z_satdy_XLAT = np.concatenate((TotlPoint_00to12Z_satdy_XLAT,extra_XLAT),axis=0)
TotlPoint_w_extra_00to12Z_satdy_CO2 = np.concatenate((TotlPoint_00to12Z_satdy_CO2,extra_CO2),axis=1)
TotlPoint_w_extra_00to12Z_satdy_CO = np.concatenate((TotlPoint_00to12Z_satdy_CO,extra_X_dict["extra_CO"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_NH3 = np.concatenate((TotlPoint_00to12Z_satdy_NH3,extra_X_dict["extra_NH3"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_NOX = np.concatenate((TotlPoint_00to12Z_satdy_NOX,extra_X_dict["extra_NOX"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM10_PRI = np.concatenate((TotlPoint_00to12Z_satdy_PM10_PRI,extra_X_dict["extra_PM10-PRI"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM25_PRI = np.concatenate((TotlPoint_00to12Z_satdy_PM25_PRI,extra_X_dict["extra_PM25-PRI"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_SO2 = np.concatenate((TotlPoint_00to12Z_satdy_SO2,extra_X_dict["extra_SO2"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_VOC = np.concatenate((TotlPoint_00to12Z_satdy_VOC,extra_X_dict["extra_VOC"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC01 = np.concatenate((TotlPoint_00to12Z_satdy_HC01,extra_X_dict["extra_HC01"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC02 = np.concatenate((TotlPoint_00to12Z_satdy_HC02,extra_X_dict["extra_HC02"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC03 = np.concatenate((TotlPoint_00to12Z_satdy_HC03,extra_X_dict["extra_HC03"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC04 = np.concatenate((TotlPoint_00to12Z_satdy_HC04,extra_X_dict["extra_HC04"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC05 = np.concatenate((TotlPoint_00to12Z_satdy_HC05,extra_X_dict["extra_HC05"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC06 = np.concatenate((TotlPoint_00to12Z_satdy_HC06,extra_X_dict["extra_HC06"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC07 = np.concatenate((TotlPoint_00to12Z_satdy_HC07,extra_X_dict["extra_HC07"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC08 = np.concatenate((TotlPoint_00to12Z_satdy_HC08,extra_X_dict["extra_HC08"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC09 = np.concatenate((TotlPoint_00to12Z_satdy_HC09,extra_X_dict["extra_HC09"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC10 = np.concatenate((TotlPoint_00to12Z_satdy_HC10,extra_X_dict["extra_HC10"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC11 = np.concatenate((TotlPoint_00to12Z_satdy_HC11,extra_X_dict["extra_HC11"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC12 = np.concatenate((TotlPoint_00to12Z_satdy_HC12,extra_X_dict["extra_HC12"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC13 = np.concatenate((TotlPoint_00to12Z_satdy_HC13,extra_X_dict["extra_HC13"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC14 = np.concatenate((TotlPoint_00to12Z_satdy_HC14,extra_X_dict["extra_HC14"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC15 = np.concatenate((TotlPoint_00to12Z_satdy_HC15,extra_X_dict["extra_HC15"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC16 = np.concatenate((TotlPoint_00to12Z_satdy_HC16,extra_X_dict["extra_HC16"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC17 = np.concatenate((TotlPoint_00to12Z_satdy_HC17,extra_X_dict["extra_HC17"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC18 = np.concatenate((TotlPoint_00to12Z_satdy_HC18,extra_X_dict["extra_HC18"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC19 = np.concatenate((TotlPoint_00to12Z_satdy_HC19,extra_X_dict["extra_HC19"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC20 = np.concatenate((TotlPoint_00to12Z_satdy_HC20,extra_X_dict["extra_HC20"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC21 = np.concatenate((TotlPoint_00to12Z_satdy_HC21,extra_X_dict["extra_HC21"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC22 = np.concatenate((TotlPoint_00to12Z_satdy_HC22,extra_X_dict["extra_HC22"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC23 = np.concatenate((TotlPoint_00to12Z_satdy_HC23,extra_X_dict["extra_HC23"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC24 = np.concatenate((TotlPoint_00to12Z_satdy_HC24,extra_X_dict["extra_HC24"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC25 = np.concatenate((TotlPoint_00to12Z_satdy_HC25,extra_X_dict["extra_HC25"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC26 = np.concatenate((TotlPoint_00to12Z_satdy_HC26,extra_X_dict["extra_HC26"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC27 = np.concatenate((TotlPoint_00to12Z_satdy_HC27,extra_X_dict["extra_HC27"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC28 = np.concatenate((TotlPoint_00to12Z_satdy_HC28,extra_X_dict["extra_HC28"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC29 = np.concatenate((TotlPoint_00to12Z_satdy_HC29,extra_X_dict["extra_HC29"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC30 = np.concatenate((TotlPoint_00to12Z_satdy_HC30,extra_X_dict["extra_HC30"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC31 = np.concatenate((TotlPoint_00to12Z_satdy_HC31,extra_X_dict["extra_HC31"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC32 = np.concatenate((TotlPoint_00to12Z_satdy_HC32,extra_X_dict["extra_HC32"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC33 = np.concatenate((TotlPoint_00to12Z_satdy_HC33,extra_X_dict["extra_HC33"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC34 = np.concatenate((TotlPoint_00to12Z_satdy_HC34,extra_X_dict["extra_HC34"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC35 = np.concatenate((TotlPoint_00to12Z_satdy_HC35,extra_X_dict["extra_HC35"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC36 = np.concatenate((TotlPoint_00to12Z_satdy_HC36,extra_X_dict["extra_HC36"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC37 = np.concatenate((TotlPoint_00to12Z_satdy_HC37,extra_X_dict["extra_HC37"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC38 = np.concatenate((TotlPoint_00to12Z_satdy_HC38,extra_X_dict["extra_HC38"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC39 = np.concatenate((TotlPoint_00to12Z_satdy_HC39,extra_X_dict["extra_HC39"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC40 = np.concatenate((TotlPoint_00to12Z_satdy_HC40,extra_X_dict["extra_HC40"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC41 = np.concatenate((TotlPoint_00to12Z_satdy_HC41,extra_X_dict["extra_HC41"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC42 = np.concatenate((TotlPoint_00to12Z_satdy_HC42,extra_X_dict["extra_HC42"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC43 = np.concatenate((TotlPoint_00to12Z_satdy_HC43,extra_X_dict["extra_HC43"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC44 = np.concatenate((TotlPoint_00to12Z_satdy_HC44,extra_X_dict["extra_HC44"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC45 = np.concatenate((TotlPoint_00to12Z_satdy_HC45,extra_X_dict["extra_HC45"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC46 = np.concatenate((TotlPoint_00to12Z_satdy_HC46,extra_X_dict["extra_HC46"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC47 = np.concatenate((TotlPoint_00to12Z_satdy_HC47,extra_X_dict["extra_HC47"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC48 = np.concatenate((TotlPoint_00to12Z_satdy_HC48,extra_X_dict["extra_HC48"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC49 = np.concatenate((TotlPoint_00to12Z_satdy_HC49,extra_X_dict["extra_HC49"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC50 = np.concatenate((TotlPoint_00to12Z_satdy_HC50,extra_X_dict["extra_HC50"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC51 = np.concatenate((TotlPoint_00to12Z_satdy_HC51,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC52 = np.concatenate((TotlPoint_00to12Z_satdy_HC52,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC53 = np.concatenate((TotlPoint_00to12Z_satdy_HC53,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC54 = np.concatenate((TotlPoint_00to12Z_satdy_HC54,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC55 = np.concatenate((TotlPoint_00to12Z_satdy_HC55,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC56 = np.concatenate((TotlPoint_00to12Z_satdy_HC56,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC57 = np.concatenate((TotlPoint_00to12Z_satdy_HC57,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC58 = np.concatenate((TotlPoint_00to12Z_satdy_HC58,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC59 = np.concatenate((TotlPoint_00to12Z_satdy_HC59,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC60 = np.concatenate((TotlPoint_00to12Z_satdy_HC60,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC61 = np.concatenate((TotlPoint_00to12Z_satdy_HC61,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC62 = np.concatenate((TotlPoint_00to12Z_satdy_HC62,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC63 = np.concatenate((TotlPoint_00to12Z_satdy_HC63,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC64 = np.concatenate((TotlPoint_00to12Z_satdy_HC64,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC65 = np.concatenate((TotlPoint_00to12Z_satdy_HC65,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC66 = np.concatenate((TotlPoint_00to12Z_satdy_HC66,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC67 = np.concatenate((TotlPoint_00to12Z_satdy_HC67,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_HC68 = np.concatenate((TotlPoint_00to12Z_satdy_HC68,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM01 = np.concatenate((TotlPoint_00to12Z_satdy_PM01,extra_X_dict["extra_PM01"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM02 = np.concatenate((TotlPoint_00to12Z_satdy_PM02,extra_X_dict["extra_PM02"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM03 = np.concatenate((TotlPoint_00to12Z_satdy_PM03,extra_X_dict["extra_PM03"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM04 = np.concatenate((TotlPoint_00to12Z_satdy_PM04,extra_X_dict["extra_PM04"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM05 = np.concatenate((TotlPoint_00to12Z_satdy_PM05,extra_X_dict["extra_PM05"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM06 = np.concatenate((TotlPoint_00to12Z_satdy_PM06,extra_X_dict["extra_PM06"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM07 = np.concatenate((TotlPoint_00to12Z_satdy_PM07,extra_X_dict["extra_PM07"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM08 = np.concatenate((TotlPoint_00to12Z_satdy_PM08,extra_X_dict["extra_PM08"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM09 = np.concatenate((TotlPoint_00to12Z_satdy_PM09,extra_X_dict["extra_PM09"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM10 = np.concatenate((TotlPoint_00to12Z_satdy_PM10,extra_X_dict["extra_PM10"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM11 = np.concatenate((TotlPoint_00to12Z_satdy_PM11,extra_X_dict["extra_PM11"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM12 = np.concatenate((TotlPoint_00to12Z_satdy_PM12,extra_X_dict["extra_PM12"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM13 = np.concatenate((TotlPoint_00to12Z_satdy_PM13,extra_X_dict["extra_PM13"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM14 = np.concatenate((TotlPoint_00to12Z_satdy_PM14,extra_X_dict["extra_PM14"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM15 = np.concatenate((TotlPoint_00to12Z_satdy_PM15,extra_X_dict["extra_PM15"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM16 = np.concatenate((TotlPoint_00to12Z_satdy_PM16,extra_X_dict["extra_PM16"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM17 = np.concatenate((TotlPoint_00to12Z_satdy_PM17,extra_X_dict["extra_PM17"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM18 = np.concatenate((TotlPoint_00to12Z_satdy_PM18,extra_X_dict["extra_PM18"]),axis=1)
TotlPoint_w_extra_00to12Z_satdy_PM19 = np.concatenate((TotlPoint_00to12Z_satdy_PM19,extra_X_dict["extra_PM19"]),axis=1)

###################################################################################################
#write total points with extra points appended
TotlPoint_w_extra_00to12Z_satdy_fn = append_dir+'/satdy/PtOG_00to12Z.nc'
TotlPoint_w_extra_00to12Z_satdy_file = Dataset(TotlPoint_w_extra_00to12Z_satdy_fn,mode='w',format='NETCDF3_64BIT')

#Creat dimensions
TotlPoint_w_extra_00to12Z_satdy_file.createDimension("ROW", nROW)
TotlPoint_w_extra_00to12Z_satdy_file.createDimension("Time", 12)
TotlPoint_w_extra_00to12Z_satdy_file.sync()

#Create variables
#float ITYPE(ROW) ;
ITYPE = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('ITYPE','f4',('ROW'),fill_value = float(0))
ITYPE[:] = TotlPoint_w_extra_00to12Z_satdy_ITYPE
varattrs=["FieldType","MemoryOrder","description","units","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['ITYPE'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['ITYPE'], varattr);
        setattr(ITYPE, varattr, varattrVal)

#float STKht(ROW) ;
STKht = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('STKht','f4',('ROW'),fill_value = float(0))
STKht[:] = TotlPoint_w_extra_00to12Z_satdy_STKht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['STKht'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['STKht'], varattr);
        setattr(STKht, varattr, varattrVal)
        
#float STKdiam(ROW) ;
STKdiam = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('STKdiam','f4',('ROW'),fill_value = float(0))
STKdiam[:] = TotlPoint_w_extra_00to12Z_satdy_STKdiam
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['STKdiam'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['STKdiam'], varattr);
        setattr(STKdiam, varattr, varattrVal)

#float STKtemp(ROW) ;
STKtemp = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('STKtemp','f4',('ROW'),fill_value = float(0))
STKtemp[:] = TotlPoint_w_extra_00to12Z_satdy_STKtemp
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['STKtemp'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['STKtemp'], varattr);
        setattr(STKtemp, varattr, varattrVal)
        
#float STKve(ROW) ;
STKve = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('STKve','f4',('ROW'),fill_value = float(0))
STKve[:] = TotlPoint_w_extra_00to12Z_satdy_STKve
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['STKve'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['STKve'], varattr);
        setattr(STKve, varattr, varattrVal)
        
#float STKflw(ROW) ;
STKflw = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('STKflw','f4',('ROW'),fill_value = float(0))
STKflw[:] = TotlPoint_w_extra_00to12Z_satdy_STKflw
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['STKflw'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['STKflw'], varattr);
        setattr(STKflw, varattr, varattrVal)
        
#float FUGht(ROW) ;
FUGht = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('FUGht','f4',('ROW'),fill_value = float(0))
FUGht[:] = TotlPoint_w_extra_00to12Z_satdy_FUGht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['FUGht'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['FUGht'], varattr);
        setattr(FUGht, varattr, varattrVal)
        
#float XLONG(ROW) ;
XLONG = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('XLONG','f4',('ROW'),fill_value = float(0))
XLONG[:] = TotlPoint_w_extra_00to12Z_satdy_XLONG
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['XLONG'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['XLONG'], varattr);
        setattr(XLONG, varattr, varattrVal)
        
#float XLAT(ROW) ;
XLAT = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('XLAT','f4',('ROW'),fill_value = float(0))
XLAT[:] = TotlPoint_w_extra_00to12Z_satdy_XLAT
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['XLAT'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['XLAT'], varattr);
        setattr(XLAT, varattr, varattrVal)

#float CO2(Time, ROW) ;
CO2 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('CO2','f4',('Time','ROW'))
CO2[:] = TotlPoint_w_extra_00to12Z_satdy_CO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['CO2'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['CO2'], varattr);
        setattr(CO2, varattr, varattrVal)
        
#float CO(Time, ROW) ;
CO = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('CO','f4',('Time','ROW'))
CO[:] = TotlPoint_w_extra_00to12Z_satdy_CO
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['CO'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['CO'], varattr);
        setattr(CO, varattr, varattrVal)
        
#float NH3(Time, ROW) ;
NH3 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('NH3','f4',('Time','ROW'))
NH3[:] = TotlPoint_w_extra_00to12Z_satdy_NH3
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['NH3'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['NH3'], varattr);
        setattr(NH3, varattr, varattrVal)
        
#float NOX(Time, ROW) ;
NOX = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('NOX','f4',('Time','ROW'))
NOX[:] = TotlPoint_w_extra_00to12Z_satdy_NOX
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['NOX'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['NOX'], varattr);
        setattr(NOX, varattr, varattrVal)
        
#float PM10-PRI(Time, ROW) ;
PM10_PRI = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM10-PRI','f4',('Time','ROW'))
PM10_PRI[:] = TotlPoint_w_extra_00to12Z_satdy_PM10_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM10-PRI'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM10-PRI'], varattr);
        setattr(PM10_PRI, varattr, varattrVal)
        
#float PM25-PRI(Time, ROW) ;
PM25_PRI = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM25-PRI','f4',('Time','ROW'))
PM25_PRI[:] = TotlPoint_w_extra_00to12Z_satdy_PM25_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM25-PRI'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM25-PRI'], varattr);
        setattr(PM25_PRI, varattr, varattrVal)
        
#float SO2(Time, ROW) ;
SO2 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('SO2','f4',('Time','ROW'))
SO2[:] = TotlPoint_w_extra_00to12Z_satdy_SO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['SO2'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['SO2'], varattr);
        setattr(SO2, varattr, varattrVal)
        
#float VOC(Time, ROW) ;
VOC = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('VOC','f4',('Time','ROW'))
VOC[:] = TotlPoint_w_extra_00to12Z_satdy_VOC
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['VOC'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['VOC'], varattr);
        setattr(VOC, varattr, varattrVal)
        
#float HC01(Time, ROW) ;
HC01 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC01','f4',('Time','ROW'))
HC01[:] = TotlPoint_w_extra_00to12Z_satdy_HC01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC01'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC01'], varattr);
        setattr(HC01, varattr, varattrVal)
        
#float HC02(Time, ROW) ;
HC02 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC02','f4',('Time','ROW'))
HC02[:] = TotlPoint_w_extra_00to12Z_satdy_HC02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC02'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC02'], varattr);
        setattr(HC02, varattr, varattrVal)
        
#float HC03(Time, ROW) ;
HC03 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC03','f4',('Time','ROW'))
HC03[:] = TotlPoint_w_extra_00to12Z_satdy_HC03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC03'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC03'], varattr);
        setattr(HC03, varattr, varattrVal)
        
#float HC04(Time, ROW) ;
HC04 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC04','f4',('Time','ROW'))
HC04[:] = TotlPoint_w_extra_00to12Z_satdy_HC04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC04'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC04'], varattr);
        setattr(HC04, varattr, varattrVal)
        
#float HC05(Time, ROW) ;
HC05 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC05','f4',('Time','ROW'))
HC05[:] = TotlPoint_w_extra_00to12Z_satdy_HC05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC05'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC05'], varattr);
        setattr(HC05, varattr, varattrVal)
        
#float HC06(Time, ROW) ;
HC06 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC06','f4',('Time','ROW'))
HC06[:] = TotlPoint_w_extra_00to12Z_satdy_HC06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC06'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC06'], varattr);
        setattr(HC06, varattr, varattrVal)
        
#float HC07(Time, ROW) ;
HC07 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC07','f4',('Time','ROW'))
HC07[:] = TotlPoint_w_extra_00to12Z_satdy_HC07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC07'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC07'], varattr);
        setattr(HC07, varattr, varattrVal)
        
#float HC08(Time, ROW) ;
HC08 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC08','f4',('Time','ROW'))
HC08[:] = TotlPoint_w_extra_00to12Z_satdy_HC08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC08'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC08'], varattr);
        setattr(HC08, varattr, varattrVal)
        
#float HC09(Time, ROW) ;
HC09 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC09','f4',('Time','ROW'))
HC09[:] = TotlPoint_w_extra_00to12Z_satdy_HC09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC09'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC09'], varattr);
        setattr(HC09, varattr, varattrVal)
        
#float HC10(Time, ROW) ;
HC10 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC10','f4',('Time','ROW'))
HC10[:] = TotlPoint_w_extra_00to12Z_satdy_HC10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC10'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC10'], varattr);
        setattr(HC10, varattr, varattrVal)

#float HC11(Time, ROW) ;
HC11 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC11','f4',('Time','ROW'))
HC11[:] = TotlPoint_w_extra_00to12Z_satdy_HC11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC11'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC11'], varattr);
        setattr(HC11, varattr, varattrVal)
        
#float HC12(Time, ROW) ;
HC12 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC12','f4',('Time','ROW'))
HC12[:] = TotlPoint_w_extra_00to12Z_satdy_HC12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC12'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC12'], varattr);
        setattr(HC12, varattr, varattrVal)
        
#float HC13(Time, ROW) ;
HC13 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC13','f4',('Time','ROW'))
HC13[:] = TotlPoint_w_extra_00to12Z_satdy_HC13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC13'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC13'], varattr);
        setattr(HC13, varattr, varattrVal)
        
#float HC14(Time, ROW) ;
HC14 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC14','f4',('Time','ROW'))
HC14[:] = TotlPoint_w_extra_00to12Z_satdy_HC14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC14'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC14'], varattr);
        setattr(HC14, varattr, varattrVal)
        
#float HC15(Time, ROW) ;
HC15 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC15','f4',('Time','ROW'))
HC15[:] = TotlPoint_w_extra_00to12Z_satdy_HC15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC15'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC15'], varattr);
        setattr(HC15, varattr, varattrVal)
        
#float HC16(Time, ROW) ;
HC16 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC16','f4',('Time','ROW'))
HC16[:] = TotlPoint_w_extra_00to12Z_satdy_HC16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC16'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC16'], varattr);
        setattr(HC16, varattr, varattrVal)
        
#float HC17(Time, ROW) ;
HC17 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC17','f4',('Time','ROW'))
HC17[:] = TotlPoint_w_extra_00to12Z_satdy_HC17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC17'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC17'], varattr);
        setattr(HC17, varattr, varattrVal)
        
#float HC18(Time, ROW) ;
HC18 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC18','f4',('Time','ROW'))
HC18[:] = TotlPoint_w_extra_00to12Z_satdy_HC18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC18'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC18'], varattr);
        setattr(HC18, varattr, varattrVal)
        
#float HC19(Time, ROW) ;
HC19 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC19','f4',('Time','ROW'))
HC19[:] = TotlPoint_w_extra_00to12Z_satdy_HC19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC19'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC19'], varattr);
        setattr(HC19, varattr, varattrVal)

#float HC20(Time, ROW) ;
HC20 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC20','f4',('Time','ROW'))
HC20[:] = TotlPoint_w_extra_00to12Z_satdy_HC20
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC20'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC20'], varattr);
        setattr(HC20, varattr, varattrVal)

#float HC21(Time, ROW) ;
HC21 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC21','f4',('Time','ROW'))
HC21[:] = TotlPoint_w_extra_00to12Z_satdy_HC21
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC21'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC21'], varattr);
        setattr(HC21, varattr, varattrVal)
        
#float HC22(Time, ROW) ;
HC22 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC22','f4',('Time','ROW'))
HC22[:] = TotlPoint_w_extra_00to12Z_satdy_HC22
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC22'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC22'], varattr);
        setattr(HC22, varattr, varattrVal)
        
#float HC23(Time, ROW) ;
HC23 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC23','f4',('Time','ROW'))
HC23[:] = TotlPoint_w_extra_00to12Z_satdy_HC23
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC23'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC23'], varattr);
        setattr(HC23, varattr, varattrVal)
        
#float HC24(Time, ROW) ;
HC24 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC24','f4',('Time','ROW'))
HC24[:] = TotlPoint_w_extra_00to12Z_satdy_HC24
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC24'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC24'], varattr);
        setattr(HC24, varattr, varattrVal)
        
#float HC25(Time, ROW) ;
HC25 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC25','f4',('Time','ROW'))
HC25[:] = TotlPoint_w_extra_00to12Z_satdy_HC25
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC25'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC25'], varattr);
        setattr(HC25, varattr, varattrVal)
        
#float HC26(Time, ROW) ;
HC26 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC26','f4',('Time','ROW'))
HC26[:] = TotlPoint_w_extra_00to12Z_satdy_HC26
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC26'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC26'], varattr);
        setattr(HC26, varattr, varattrVal)
        
#float HC27(Time, ROW) ;
HC27 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC27','f4',('Time','ROW'))
HC27[:] = TotlPoint_w_extra_00to12Z_satdy_HC27
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC27'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC27'], varattr);
        setattr(HC27, varattr, varattrVal)
        
#float HC28(Time, ROW) ;
HC28 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC28','f4',('Time','ROW'))
HC28[:] = TotlPoint_w_extra_00to12Z_satdy_HC28
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC28'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC28'], varattr);
        setattr(HC28, varattr, varattrVal)
        
#float HC29(Time, ROW) ;
HC29 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC29','f4',('Time','ROW'))
HC29[:] = TotlPoint_w_extra_00to12Z_satdy_HC29
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC29'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC29'], varattr);
        setattr(HC29, varattr, varattrVal)

#float HC30(Time, ROW) ;
HC30 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC30','f4',('Time','ROW'))
HC30[:] = TotlPoint_w_extra_00to12Z_satdy_HC30
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC30'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC30'], varattr);
        setattr(HC30, varattr, varattrVal)

#float HC31(Time, ROW) ;
HC31 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC31','f4',('Time','ROW'))
HC31[:] = TotlPoint_w_extra_00to12Z_satdy_HC31
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC31'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC31'], varattr);
        setattr(HC31, varattr, varattrVal)
        
#float HC32(Time, ROW) ;
HC32 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC32','f4',('Time','ROW'))
HC32[:] = TotlPoint_w_extra_00to12Z_satdy_HC32
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC32'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC32'], varattr);
        setattr(HC32, varattr, varattrVal)
        
#float HC33(Time, ROW) ;
HC33 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC33','f4',('Time','ROW'))
HC33[:] = TotlPoint_w_extra_00to12Z_satdy_HC33
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC33'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC33'], varattr);
        setattr(HC33, varattr, varattrVal)
        
#float HC34(Time, ROW) ;
HC34 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC34','f4',('Time','ROW'))
HC34[:] = TotlPoint_w_extra_00to12Z_satdy_HC34
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC34'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC34'], varattr);
        setattr(HC34, varattr, varattrVal)
        
#float HC35(Time, ROW) ;
HC35 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC35','f4',('Time','ROW'))
HC35[:] = TotlPoint_w_extra_00to12Z_satdy_HC35
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC35'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC35'], varattr);
        setattr(HC35, varattr, varattrVal)
        
#float HC36(Time, ROW) ;
HC36 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC36','f4',('Time','ROW'))
HC36[:] = TotlPoint_w_extra_00to12Z_satdy_HC36
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC36'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC36'], varattr);
        setattr(HC36, varattr, varattrVal)
        
#float HC37(Time, ROW) ;
HC37 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC37','f4',('Time','ROW'))
HC37[:] = TotlPoint_w_extra_00to12Z_satdy_HC37
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC37'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC37'], varattr);
        setattr(HC37, varattr, varattrVal)
        
#float HC38(Time, ROW) ;
HC38 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC38','f4',('Time','ROW'))
HC38[:] = TotlPoint_w_extra_00to12Z_satdy_HC38
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC38'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC38'], varattr);
        setattr(HC38, varattr, varattrVal)
        
#float HC39(Time, ROW) ;
HC39 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC39','f4',('Time','ROW'))
HC39[:] = TotlPoint_w_extra_00to12Z_satdy_HC39
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC39'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC39'], varattr);
        setattr(HC39, varattr, varattrVal)
        
#float HC40(Time, ROW) ;
HC40 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC40','f4',('Time','ROW'))
HC40[:] = TotlPoint_w_extra_00to12Z_satdy_HC40
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC40'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC40'], varattr);
        setattr(HC40, varattr, varattrVal)

#float HC41(Time, ROW) ;
HC41 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC41','f4',('Time','ROW'))
HC41[:] = TotlPoint_w_extra_00to12Z_satdy_HC41
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC41'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC41'], varattr);
        setattr(HC41, varattr, varattrVal)
        
#float HC42(Time, ROW) ;
HC42 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC42','f4',('Time','ROW'))
HC42[:] = TotlPoint_w_extra_00to12Z_satdy_HC42
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC42'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC42'], varattr);
        setattr(HC42, varattr, varattrVal)
        
#float HC43(Time, ROW) ;
HC43 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC43','f4',('Time','ROW'))
HC43[:] = TotlPoint_w_extra_00to12Z_satdy_HC43
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC43'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC43'], varattr);
        setattr(HC43, varattr, varattrVal)
        
#float HC44(Time, ROW) ;
HC44 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC44','f4',('Time','ROW'))
HC44[:] = TotlPoint_w_extra_00to12Z_satdy_HC44
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC44'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC44'], varattr);
        setattr(HC44, varattr, varattrVal)
        
#float HC45(Time, ROW) ;
HC45 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC45','f4',('Time','ROW'))
HC45[:] = TotlPoint_w_extra_00to12Z_satdy_HC45
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC45'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC45'], varattr);
        setattr(HC45, varattr, varattrVal)
        
#float HC46(Time, ROW) ;
HC46 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC46','f4',('Time','ROW'))
HC46[:] = TotlPoint_w_extra_00to12Z_satdy_HC46
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC46'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC46'], varattr);
        setattr(HC46, varattr, varattrVal)
        
#float HC47(Time, ROW) ;
HC47 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC47','f4',('Time','ROW'))
HC47[:] = TotlPoint_w_extra_00to12Z_satdy_HC47
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC47'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC47'], varattr);
        setattr(HC47, varattr, varattrVal)
        
#float HC48(Time, ROW) ;
HC48 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC48','f4',('Time','ROW'))
HC48[:] = TotlPoint_w_extra_00to12Z_satdy_HC48
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC48'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC48'], varattr);
        setattr(HC48, varattr, varattrVal)
        
#float HC49(Time, ROW) ;
HC49 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC49','f4',('Time','ROW'))
HC49[:] = TotlPoint_w_extra_00to12Z_satdy_HC49
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC49'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC49'], varattr);
        setattr(HC49, varattr, varattrVal)
        
#float HC50(Time, ROW) ;
HC50 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC50','f4',('Time','ROW'))
HC50[:] = TotlPoint_w_extra_00to12Z_satdy_HC50
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC50'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC50'], varattr);
        setattr(HC50, varattr, varattrVal)

#float HC51(Time, ROW) ;
HC51 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC51','f4',('Time','ROW'))
HC51[:] = TotlPoint_w_extra_00to12Z_satdy_HC51
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC51'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC51'], varattr);
        setattr(HC51, varattr, varattrVal)
        
#float HC52(Time, ROW) ;
HC52 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC52','f4',('Time','ROW'))
HC52[:] = TotlPoint_w_extra_00to12Z_satdy_HC52
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC52'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC52'], varattr);
        setattr(HC52, varattr, varattrVal)
        
#float HC53(Time, ROW) ;
HC53 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC53','f4',('Time','ROW'))
HC53[:] = TotlPoint_w_extra_00to12Z_satdy_HC53
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC53'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC53'], varattr);
        setattr(HC53, varattr, varattrVal)
        
#float HC54(Time, ROW) ;
HC54 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC54','f4',('Time','ROW'))
HC54[:] = TotlPoint_w_extra_00to12Z_satdy_HC54
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC54'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC54'], varattr);
        setattr(HC54, varattr, varattrVal)
        
#float HC55(Time, ROW) ;
HC55 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC55','f4',('Time','ROW'))
HC55[:] = TotlPoint_w_extra_00to12Z_satdy_HC55
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC55'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC55'], varattr);
        setattr(HC55, varattr, varattrVal)
        
#float HC56(Time, ROW) ;
HC56 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC56','f4',('Time','ROW'))
HC56[:] = TotlPoint_w_extra_00to12Z_satdy_HC56
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC56'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC56'], varattr);
        setattr(HC56, varattr, varattrVal)
        
#float HC57(Time, ROW) ;
HC57 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC57','f4',('Time','ROW'))
HC57[:] = TotlPoint_w_extra_00to12Z_satdy_HC57
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC57'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC57'], varattr);
        setattr(HC57, varattr, varattrVal)
        
#float HC58(Time, ROW) ;
HC58 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC58','f4',('Time','ROW'))
HC58[:] = TotlPoint_w_extra_00to12Z_satdy_HC58
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC58'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC58'], varattr);
        setattr(HC58, varattr, varattrVal)
        
#float HC59(Time, ROW) ;
HC59 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC59','f4',('Time','ROW'))
HC59[:] = TotlPoint_w_extra_00to12Z_satdy_HC59
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC59'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC59'], varattr);
        setattr(HC59, varattr, varattrVal)

#float HC60(Time, ROW) ;
HC60 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC60','f4',('Time','ROW'))
HC60[:] = TotlPoint_w_extra_00to12Z_satdy_HC60
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC60'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC60'], varattr);
        setattr(HC60, varattr, varattrVal)

#float HC61(Time, ROW) ;
HC61 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC61','f4',('Time','ROW'))
HC61[:] = TotlPoint_w_extra_00to12Z_satdy_HC61
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC61'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC61'], varattr);
        setattr(HC61, varattr, varattrVal)
        
#float HC62(Time, ROW) ;
HC62 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC62','f4',('Time','ROW'))
HC62[:] = TotlPoint_w_extra_00to12Z_satdy_HC62
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC62'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC62'], varattr);
        setattr(HC62, varattr, varattrVal)
        
#float HC63(Time, ROW) ;
HC63 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC63','f4',('Time','ROW'))
HC63[:] = TotlPoint_w_extra_00to12Z_satdy_HC63
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC63'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC63'], varattr);
        setattr(HC63, varattr, varattrVal)
        
#float HC64(Time, ROW) ;
HC64 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC64','f4',('Time','ROW'))
HC64[:] = TotlPoint_w_extra_00to12Z_satdy_HC64
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC64'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC64'], varattr);
        setattr(HC64, varattr, varattrVal)
        
#float HC65(Time, ROW) ;
HC65 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC65','f4',('Time','ROW'))
HC65[:] = TotlPoint_w_extra_00to12Z_satdy_HC65
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC65'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC65'], varattr);
        setattr(HC65, varattr, varattrVal)
        
#float HC66(Time, ROW) ;
HC66 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC66','f4',('Time','ROW'))
HC66[:] = TotlPoint_w_extra_00to12Z_satdy_HC66
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC66'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC66'], varattr);
        setattr(HC66, varattr, varattrVal)
        
#float HC67(Time, ROW) ;
HC67 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC67','f4',('Time','ROW'))
HC67[:] = TotlPoint_w_extra_00to12Z_satdy_HC67
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC67'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC67'], varattr);
        setattr(HC67, varattr, varattrVal)
        
#float HC68(Time, ROW) ;
HC68 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('HC68','f4',('Time','ROW'))
HC68[:] = TotlPoint_w_extra_00to12Z_satdy_HC68
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['HC68'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['HC68'], varattr);
        setattr(HC68, varattr, varattrVal)

#float PM01(Time, ROW) ;
PM01 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM01','f4',('Time','ROW'))
PM01[:] = TotlPoint_w_extra_00to12Z_satdy_PM01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM01'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM01'], varattr);
        setattr(PM01, varattr, varattrVal)

#float PM02(Time, ROW) ;
PM02 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM02','f4',('Time','ROW'))
PM02[:] = TotlPoint_w_extra_00to12Z_satdy_PM02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM02'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM02'], varattr);
        setattr(PM02, varattr, varattrVal)
        
#float PM03(Time, ROW) ;
PM03 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM03','f4',('Time','ROW'))
PM03[:] = TotlPoint_w_extra_00to12Z_satdy_PM03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM03'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM03'], varattr);
        setattr(PM03, varattr, varattrVal)

#float PM04(Time, ROW) ;
PM04 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM04','f4',('Time','ROW'))
PM04[:] = TotlPoint_w_extra_00to12Z_satdy_PM04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM04'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM04'], varattr);
        setattr(PM04, varattr, varattrVal)
        
#float PM05(Time, ROW) ;
PM05 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM05','f4',('Time','ROW'))
PM05[:] = TotlPoint_w_extra_00to12Z_satdy_PM05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM05'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM05'], varattr);
        setattr(PM05, varattr, varattrVal)

#float PM06(Time, ROW) ;
PM06 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM06','f4',('Time','ROW'))
PM06[:] = TotlPoint_w_extra_00to12Z_satdy_PM06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM06'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM06'], varattr);
        setattr(PM06, varattr, varattrVal)
        
#float PM07(Time, ROW) ;
PM07 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM07','f4',('Time','ROW'))
PM07[:] = TotlPoint_w_extra_00to12Z_satdy_PM07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM07'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM07'], varattr);
        setattr(PM07, varattr, varattrVal)
        
#float PM08(Time, ROW) ;
PM08 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM08','f4',('Time','ROW'))
PM08[:] = TotlPoint_w_extra_00to12Z_satdy_PM08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM08'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM08'], varattr);
        setattr(PM08, varattr, varattrVal)
        
#float PM09(Time, ROW) ;
PM09 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM09','f4',('Time','ROW'))
PM09[:] = TotlPoint_w_extra_00to12Z_satdy_PM09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM09'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM09'], varattr);
        setattr(PM09, varattr, varattrVal)
        
#float PM10(Time, ROW) ;
PM10 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM10','f4',('Time','ROW'))
PM10[:] = TotlPoint_w_extra_00to12Z_satdy_PM10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM10'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM10'], varattr);
        setattr(PM10, varattr, varattrVal)
        
#float PM11(Time, ROW) ;
PM11 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM11','f4',('Time','ROW'))
PM11[:] = TotlPoint_w_extra_00to12Z_satdy_PM11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM11'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM11'], varattr);
        setattr(PM11, varattr, varattrVal)

#float PM12(Time, ROW) ;
PM12 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM12','f4',('Time','ROW'))
PM12[:] = TotlPoint_w_extra_00to12Z_satdy_PM12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM12'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM12'], varattr);
        setattr(PM12, varattr, varattrVal)
        
#float PM13(Time, ROW) ;
PM13 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM13','f4',('Time','ROW'))
PM13[:] = TotlPoint_w_extra_00to12Z_satdy_PM13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM13'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM13'], varattr);
        setattr(PM13, varattr, varattrVal)

#float PM14(Time, ROW) ;
PM14 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM14','f4',('Time','ROW'))
PM14[:] = TotlPoint_w_extra_00to12Z_satdy_PM14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM14'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM14'], varattr);
        setattr(PM14, varattr, varattrVal)
        
#float PM15(Time, ROW) ;
PM15 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM15','f4',('Time','ROW'))
PM15[:] = TotlPoint_w_extra_00to12Z_satdy_PM15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM15'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM15'], varattr);
        setattr(PM15, varattr, varattrVal)

#float PM16(Time, ROW) ;
PM16 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM16','f4',('Time','ROW'))
PM16[:] = TotlPoint_w_extra_00to12Z_satdy_PM16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM16'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM16'], varattr);
        setattr(PM16, varattr, varattrVal)
        
#float PM17(Time, ROW) ;
PM17 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM17','f4',('Time','ROW'))
PM17[:] = TotlPoint_w_extra_00to12Z_satdy_PM17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM17'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM17'], varattr);
        setattr(PM17, varattr, varattrVal)
        
#float PM18(Time, ROW) ;
PM18 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM18','f4',('Time','ROW'))
PM18[:] = TotlPoint_w_extra_00to12Z_satdy_PM18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM18'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM18'], varattr);
        setattr(PM18, varattr, varattrVal)
        
#float PM19(Time, ROW) ;
PM19 = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('PM19','f4',('Time','ROW'))
PM19[:] = TotlPoint_w_extra_00to12Z_satdy_PM19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_satdy_file.variables['PM19'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file.variables['PM19'], varattr);
        setattr(PM19, varattr, varattrVal)

#char Times(Time) ;
Times = TotlPoint_w_extra_00to12Z_satdy_file.createVariable('Times','S1',('Time'))
Times[:] = TotlPoint_00to12Z_satdy_Times

#copy global attributes from TotlPoint_00to12Z_satdy_file
for varattr in TotlPoint_00to12Z_satdy_file.ncattrs():
    if hasattr(TotlPoint_00to12Z_satdy_file, varattr):
        varattrVal = getattr(TotlPoint_00to12Z_satdy_file, varattr);
        setattr(TotlPoint_w_extra_00to12Z_satdy_file, varattr, varattrVal)

TotlPoint_w_extra_00to12Z_satdy_file.close()


# In[19]:


#append extra points to original point file that is input to wrfchemi assembly program

###################################################################################################
#satdy, 12to24Z

###################################################################################################
#read original variables
TotlPoint_12to24Z_satdy_fn = base_dir+'/satdy/PtOG_12to24Z.nc'
TotlPoint_12to24Z_satdy_file = Dataset(TotlPoint_12to24Z_satdy_fn,mode='r',open=True)
TotlPoint_12to24Z_satdy_ITYPE = TotlPoint_12to24Z_satdy_file.variables['ITYPE'][:]
TotlPoint_12to24Z_satdy_STKht = TotlPoint_12to24Z_satdy_file.variables['STKht'][:]
TotlPoint_12to24Z_satdy_STKdiam = TotlPoint_12to24Z_satdy_file.variables['STKdiam'][:]
TotlPoint_12to24Z_satdy_STKtemp = TotlPoint_12to24Z_satdy_file.variables['STKtemp'][:]
TotlPoint_12to24Z_satdy_STKve = TotlPoint_12to24Z_satdy_file.variables['STKve'][:]
TotlPoint_12to24Z_satdy_STKflw = TotlPoint_12to24Z_satdy_file.variables['STKflw'][:]
TotlPoint_12to24Z_satdy_FUGht = TotlPoint_12to24Z_satdy_file.variables['FUGht'][:]
TotlPoint_12to24Z_satdy_XLONG = TotlPoint_12to24Z_satdy_file.variables['XLONG'][:]
TotlPoint_12to24Z_satdy_XLAT = TotlPoint_12to24Z_satdy_file.variables['XLAT'][:]
TotlPoint_12to24Z_satdy_CO2 = TotlPoint_12to24Z_satdy_file.variables['CO2'][:][:] 
TotlPoint_12to24Z_satdy_CO = TotlPoint_12to24Z_satdy_file.variables['CO'][:][:] 
TotlPoint_12to24Z_satdy_NH3 = TotlPoint_12to24Z_satdy_file.variables['NH3'][:][:] 
TotlPoint_12to24Z_satdy_NOX = TotlPoint_12to24Z_satdy_file.variables['NOX'][:][:] 
TotlPoint_12to24Z_satdy_PM10_PRI = TotlPoint_12to24Z_satdy_file.variables['PM10-PRI'][:][:] 
TotlPoint_12to24Z_satdy_PM25_PRI = TotlPoint_12to24Z_satdy_file.variables['PM25-PRI'][:][:] 
TotlPoint_12to24Z_satdy_SO2 = TotlPoint_12to24Z_satdy_file.variables['SO2'][:][:] 
TotlPoint_12to24Z_satdy_VOC = TotlPoint_12to24Z_satdy_file.variables['VOC'][:][:] 
TotlPoint_12to24Z_satdy_HC01 = TotlPoint_12to24Z_satdy_file.variables['HC01'][:][:] 
TotlPoint_12to24Z_satdy_HC02 = TotlPoint_12to24Z_satdy_file.variables['HC02'][:][:] 
TotlPoint_12to24Z_satdy_HC03 = TotlPoint_12to24Z_satdy_file.variables['HC03'][:][:] 
TotlPoint_12to24Z_satdy_HC04 = TotlPoint_12to24Z_satdy_file.variables['HC04'][:][:] 
TotlPoint_12to24Z_satdy_HC05 = TotlPoint_12to24Z_satdy_file.variables['HC05'][:][:] 
TotlPoint_12to24Z_satdy_HC06 = TotlPoint_12to24Z_satdy_file.variables['HC06'][:][:] 
TotlPoint_12to24Z_satdy_HC07 = TotlPoint_12to24Z_satdy_file.variables['HC07'][:][:] 
TotlPoint_12to24Z_satdy_HC08 = TotlPoint_12to24Z_satdy_file.variables['HC08'][:][:] 
TotlPoint_12to24Z_satdy_HC09 = TotlPoint_12to24Z_satdy_file.variables['HC09'][:][:] 
TotlPoint_12to24Z_satdy_HC10 = TotlPoint_12to24Z_satdy_file.variables['HC10'][:][:] 
TotlPoint_12to24Z_satdy_HC11 = TotlPoint_12to24Z_satdy_file.variables['HC11'][:][:] 
TotlPoint_12to24Z_satdy_HC12 = TotlPoint_12to24Z_satdy_file.variables['HC12'][:][:] 
TotlPoint_12to24Z_satdy_HC13 = TotlPoint_12to24Z_satdy_file.variables['HC13'][:][:] 
TotlPoint_12to24Z_satdy_HC14 = TotlPoint_12to24Z_satdy_file.variables['HC14'][:][:] 
TotlPoint_12to24Z_satdy_HC15 = TotlPoint_12to24Z_satdy_file.variables['HC15'][:][:] 
TotlPoint_12to24Z_satdy_HC16 = TotlPoint_12to24Z_satdy_file.variables['HC16'][:][:] 
TotlPoint_12to24Z_satdy_HC17 = TotlPoint_12to24Z_satdy_file.variables['HC17'][:][:] 
TotlPoint_12to24Z_satdy_HC18 = TotlPoint_12to24Z_satdy_file.variables['HC18'][:][:] 
TotlPoint_12to24Z_satdy_HC19 = TotlPoint_12to24Z_satdy_file.variables['HC19'][:][:] 
TotlPoint_12to24Z_satdy_HC20 = TotlPoint_12to24Z_satdy_file.variables['HC20'][:][:] 
TotlPoint_12to24Z_satdy_HC21 = TotlPoint_12to24Z_satdy_file.variables['HC21'][:][:] 
TotlPoint_12to24Z_satdy_HC22 = TotlPoint_12to24Z_satdy_file.variables['HC22'][:][:] 
TotlPoint_12to24Z_satdy_HC23 = TotlPoint_12to24Z_satdy_file.variables['HC23'][:][:] 
TotlPoint_12to24Z_satdy_HC24 = TotlPoint_12to24Z_satdy_file.variables['HC24'][:][:] 
TotlPoint_12to24Z_satdy_HC25 = TotlPoint_12to24Z_satdy_file.variables['HC25'][:][:] 
TotlPoint_12to24Z_satdy_HC26 = TotlPoint_12to24Z_satdy_file.variables['HC26'][:][:] 
TotlPoint_12to24Z_satdy_HC27 = TotlPoint_12to24Z_satdy_file.variables['HC27'][:][:] 
TotlPoint_12to24Z_satdy_HC28 = TotlPoint_12to24Z_satdy_file.variables['HC28'][:][:] 
TotlPoint_12to24Z_satdy_HC29 = TotlPoint_12to24Z_satdy_file.variables['HC29'][:][:] 
TotlPoint_12to24Z_satdy_HC30 = TotlPoint_12to24Z_satdy_file.variables['HC30'][:][:] 
TotlPoint_12to24Z_satdy_HC31 = TotlPoint_12to24Z_satdy_file.variables['HC31'][:][:] 
TotlPoint_12to24Z_satdy_HC32 = TotlPoint_12to24Z_satdy_file.variables['HC32'][:][:] 
TotlPoint_12to24Z_satdy_HC33 = TotlPoint_12to24Z_satdy_file.variables['HC33'][:][:] 
TotlPoint_12to24Z_satdy_HC34 = TotlPoint_12to24Z_satdy_file.variables['HC34'][:][:] 
TotlPoint_12to24Z_satdy_HC35 = TotlPoint_12to24Z_satdy_file.variables['HC35'][:][:] 
TotlPoint_12to24Z_satdy_HC36 = TotlPoint_12to24Z_satdy_file.variables['HC36'][:][:] 
TotlPoint_12to24Z_satdy_HC37 = TotlPoint_12to24Z_satdy_file.variables['HC37'][:][:] 
TotlPoint_12to24Z_satdy_HC38 = TotlPoint_12to24Z_satdy_file.variables['HC38'][:][:] 
TotlPoint_12to24Z_satdy_HC39 = TotlPoint_12to24Z_satdy_file.variables['HC39'][:][:] 
TotlPoint_12to24Z_satdy_HC40 = TotlPoint_12to24Z_satdy_file.variables['HC40'][:][:] 
TotlPoint_12to24Z_satdy_HC41 = TotlPoint_12to24Z_satdy_file.variables['HC41'][:][:] 
TotlPoint_12to24Z_satdy_HC42 = TotlPoint_12to24Z_satdy_file.variables['HC42'][:][:] 
TotlPoint_12to24Z_satdy_HC43 = TotlPoint_12to24Z_satdy_file.variables['HC43'][:][:] 
TotlPoint_12to24Z_satdy_HC44 = TotlPoint_12to24Z_satdy_file.variables['HC44'][:][:] 
TotlPoint_12to24Z_satdy_HC45 = TotlPoint_12to24Z_satdy_file.variables['HC45'][:][:] 
TotlPoint_12to24Z_satdy_HC46 = TotlPoint_12to24Z_satdy_file.variables['HC46'][:][:] 
TotlPoint_12to24Z_satdy_HC47 = TotlPoint_12to24Z_satdy_file.variables['HC47'][:][:] 
TotlPoint_12to24Z_satdy_HC48 = TotlPoint_12to24Z_satdy_file.variables['HC48'][:][:] 
TotlPoint_12to24Z_satdy_HC49 = TotlPoint_12to24Z_satdy_file.variables['HC49'][:][:] 
TotlPoint_12to24Z_satdy_HC50 = TotlPoint_12to24Z_satdy_file.variables['HC50'][:][:] 
TotlPoint_12to24Z_satdy_HC51 = TotlPoint_12to24Z_satdy_file.variables['HC51'][:][:] 
TotlPoint_12to24Z_satdy_HC52 = TotlPoint_12to24Z_satdy_file.variables['HC52'][:][:] 
TotlPoint_12to24Z_satdy_HC53 = TotlPoint_12to24Z_satdy_file.variables['HC53'][:][:] 
TotlPoint_12to24Z_satdy_HC54 = TotlPoint_12to24Z_satdy_file.variables['HC54'][:][:] 
TotlPoint_12to24Z_satdy_HC55 = TotlPoint_12to24Z_satdy_file.variables['HC55'][:][:] 
TotlPoint_12to24Z_satdy_HC56 = TotlPoint_12to24Z_satdy_file.variables['HC56'][:][:] 
TotlPoint_12to24Z_satdy_HC57 = TotlPoint_12to24Z_satdy_file.variables['HC57'][:][:] 
TotlPoint_12to24Z_satdy_HC58 = TotlPoint_12to24Z_satdy_file.variables['HC58'][:][:] 
TotlPoint_12to24Z_satdy_HC59 = TotlPoint_12to24Z_satdy_file.variables['HC59'][:][:] 
TotlPoint_12to24Z_satdy_HC60 = TotlPoint_12to24Z_satdy_file.variables['HC60'][:][:] 
TotlPoint_12to24Z_satdy_HC61 = TotlPoint_12to24Z_satdy_file.variables['HC61'][:][:] 
TotlPoint_12to24Z_satdy_HC62 = TotlPoint_12to24Z_satdy_file.variables['HC62'][:][:] 
TotlPoint_12to24Z_satdy_HC63 = TotlPoint_12to24Z_satdy_file.variables['HC63'][:][:] 
TotlPoint_12to24Z_satdy_HC64 = TotlPoint_12to24Z_satdy_file.variables['HC64'][:][:] 
TotlPoint_12to24Z_satdy_HC65 = TotlPoint_12to24Z_satdy_file.variables['HC65'][:][:] 
TotlPoint_12to24Z_satdy_HC66 = TotlPoint_12to24Z_satdy_file.variables['HC66'][:][:] 
TotlPoint_12to24Z_satdy_HC67 = TotlPoint_12to24Z_satdy_file.variables['HC67'][:][:] 
TotlPoint_12to24Z_satdy_HC68 = TotlPoint_12to24Z_satdy_file.variables['HC68'][:][:] 
TotlPoint_12to24Z_satdy_PM01 = TotlPoint_12to24Z_satdy_file.variables['PM01'][:][:] 
TotlPoint_12to24Z_satdy_PM02 = TotlPoint_12to24Z_satdy_file.variables['PM02'][:][:] 
TotlPoint_12to24Z_satdy_PM03 = TotlPoint_12to24Z_satdy_file.variables['PM03'][:][:] 
TotlPoint_12to24Z_satdy_PM04 = TotlPoint_12to24Z_satdy_file.variables['PM04'][:][:] 
TotlPoint_12to24Z_satdy_PM05 = TotlPoint_12to24Z_satdy_file.variables['PM05'][:][:] 
TotlPoint_12to24Z_satdy_PM06 = TotlPoint_12to24Z_satdy_file.variables['PM06'][:][:] 
TotlPoint_12to24Z_satdy_PM07 = TotlPoint_12to24Z_satdy_file.variables['PM07'][:][:] 
TotlPoint_12to24Z_satdy_PM08 = TotlPoint_12to24Z_satdy_file.variables['PM08'][:][:] 
TotlPoint_12to24Z_satdy_PM09 = TotlPoint_12to24Z_satdy_file.variables['PM09'][:][:] 
TotlPoint_12to24Z_satdy_PM10 = TotlPoint_12to24Z_satdy_file.variables['PM10'][:][:] 
TotlPoint_12to24Z_satdy_PM11 = TotlPoint_12to24Z_satdy_file.variables['PM11'][:][:] 
TotlPoint_12to24Z_satdy_PM12 = TotlPoint_12to24Z_satdy_file.variables['PM12'][:][:] 
TotlPoint_12to24Z_satdy_PM13 = TotlPoint_12to24Z_satdy_file.variables['PM13'][:][:] 
TotlPoint_12to24Z_satdy_PM14 = TotlPoint_12to24Z_satdy_file.variables['PM14'][:][:] 
TotlPoint_12to24Z_satdy_PM15 = TotlPoint_12to24Z_satdy_file.variables['PM15'][:][:] 
TotlPoint_12to24Z_satdy_PM16 = TotlPoint_12to24Z_satdy_file.variables['PM16'][:][:] 
TotlPoint_12to24Z_satdy_PM17 = TotlPoint_12to24Z_satdy_file.variables['PM17'][:][:] 
TotlPoint_12to24Z_satdy_PM18 = TotlPoint_12to24Z_satdy_file.variables['PM18'][:][:] 
TotlPoint_12to24Z_satdy_PM19 = TotlPoint_12to24Z_satdy_file.variables['PM19'][:][:] 
TotlPoint_12to24Z_satdy_Times = TotlPoint_12to24Z_satdy_file.variables['Times'][:]

###################################################################################################
#get total ROW
nROW_org, = TotlPoint_12to24Z_satdy_ITYPE.shape
nROW_extra_OG = len(LON_ng_proc)
nROW_extra = nROW_extra_OG
nROW = nROW_org + nROW_extra
print("nROW_org",nROW_org)
print("nROW_extra",nROW_extra)
print("nROW",nROW)

###################################################################################################
#Organize extra_data
extra_ITYPE_OG = np.array(ERPTYPE_ng_proc)
extra_ITYPE = extra_ITYPE_OG

extra_STKht_OG = np.array(STKHGT_ng_proc)
extra_STKht = extra_STKht_OG

extra_STKdiam_OG = np.array(STKDIAM_ng_proc)
extra_STKdiam = extra_STKdiam_OG

extra_STKtemp_OG = np.array(STKTEMP_ng_proc)
extra_STKtemp = extra_STKtemp_OG

extra_STKve_OG = np.array(STKVEL_ng_proc)
extra_STKve = extra_STKve_OG

extra_STKflw_OG = np.array(STKFLOW_ng_proc)
extra_STKflw = extra_STKflw_OG

extra_FUGht = np.empty(nROW_extra) #FUGht set as empty

extra_XLONG_OG = np.array(LON_ng_proc)
extra_XLONG = extra_XLONG_OG

extra_XLAT_OG = np.array(LAT_ng_proc)
extra_XLAT = extra_XLAT_OG

extra_STATE_OG = STATE_ng_proc

###################################################################################################
#CO2

##################################################################################
extra_CO2_FCPE_ng_proc = HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_satdy[12:24,:]
extra_CO2_OG = extra_CO2_FCPE_ng_proc

##################################################################################
extra_CO2 =extra_CO2_OG

###################################################################################################
#CH4 from OG can use GHGRP numbers

##################################################################################
extra_CH4_FCPE_ng_proc = HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_satdy[12:24,:]
extra_CH4_OG = extra_CH4_FCPE_ng_proc

###################################################################################################

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#AQ: using state-level emission ratios to CO2

###################################################################################################
#grab emission ratios from the summary ratio arrays
###################################################################################################

#based on OG process type, species, and state location
extra_X_OG = np.empty([12,nROW_extra_OG,len(species_vector)])
print("extra_X_OG.shape", extra_X_OG.shape)

for pt in range(0,nROW_extra_OG):
    proc_index = 0
    #print("proc_index",proc_index)

    state_abb_cur = extra_STATE_OG[pt]
    
    if state_abb_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_abb_cur)
        #print("state_index",state_index)
    
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * proc_spec_state_emisXdCO2_OG[proc_index,spec_index,state_index]

    else:
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * statistics.mean(proc_spec_state_emisXdCO2_OG[proc_index,spec_index,:])

extra_X_OG_dict = {}
for spec_cur in species_vector:
    #for CH4 use GHGRP numbers
    if spec_cur == 'HC01':
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_CH4_OG
    else:
        spec_index = species_vector.index(spec_cur)
        #print("spec_index",spec_index)
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_X_OG[:,:,spec_index]

###################################################################################################
#stack OG AQ species
extra_X_dict = {}
for spec_cur in species_vector:
    extra_Xi_OG = extra_X_OG_dict["extra_{0}_OG".format(spec_cur)]
    extra_X = extra_Xi_OG
    extra_X_dict["extra_{0}".format(spec_cur)] = extra_X

###################################################################################################
extra_other_spec = np.zeros([12,nROW_extra])

###################################################################################################
#append extra points to original data
TotlPoint_w_extra_12to24Z_satdy_ITYPE = np.concatenate((TotlPoint_12to24Z_satdy_ITYPE,extra_ITYPE),axis=0)
TotlPoint_w_extra_12to24Z_satdy_STKht = np.concatenate((TotlPoint_12to24Z_satdy_STKht,extra_STKht),axis=0)
TotlPoint_w_extra_12to24Z_satdy_STKdiam = np.concatenate((TotlPoint_12to24Z_satdy_STKdiam,extra_STKdiam),axis=0)
TotlPoint_w_extra_12to24Z_satdy_STKtemp = np.concatenate((TotlPoint_12to24Z_satdy_STKtemp,extra_STKtemp),axis=0)
TotlPoint_w_extra_12to24Z_satdy_STKve = np.concatenate((TotlPoint_12to24Z_satdy_STKve,extra_STKve),axis=0)
TotlPoint_w_extra_12to24Z_satdy_STKflw = np.concatenate((TotlPoint_12to24Z_satdy_STKflw,extra_STKflw),axis=0)
TotlPoint_w_extra_12to24Z_satdy_FUGht = np.concatenate((TotlPoint_12to24Z_satdy_FUGht,extra_FUGht),axis=0)
TotlPoint_w_extra_12to24Z_satdy_XLONG = np.concatenate((TotlPoint_12to24Z_satdy_XLONG,extra_XLONG),axis=0)
TotlPoint_w_extra_12to24Z_satdy_XLAT = np.concatenate((TotlPoint_12to24Z_satdy_XLAT,extra_XLAT),axis=0)
TotlPoint_w_extra_12to24Z_satdy_CO2 = np.concatenate((TotlPoint_12to24Z_satdy_CO2,extra_CO2),axis=1)
TotlPoint_w_extra_12to24Z_satdy_CO = np.concatenate((TotlPoint_12to24Z_satdy_CO,extra_X_dict["extra_CO"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_NH3 = np.concatenate((TotlPoint_12to24Z_satdy_NH3,extra_X_dict["extra_NH3"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_NOX = np.concatenate((TotlPoint_12to24Z_satdy_NOX,extra_X_dict["extra_NOX"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM10_PRI = np.concatenate((TotlPoint_12to24Z_satdy_PM10_PRI,extra_X_dict["extra_PM10-PRI"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM25_PRI = np.concatenate((TotlPoint_12to24Z_satdy_PM25_PRI,extra_X_dict["extra_PM25-PRI"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_SO2 = np.concatenate((TotlPoint_12to24Z_satdy_SO2,extra_X_dict["extra_SO2"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_VOC = np.concatenate((TotlPoint_12to24Z_satdy_VOC,extra_X_dict["extra_VOC"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC01 = np.concatenate((TotlPoint_12to24Z_satdy_HC01,extra_X_dict["extra_HC01"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC02 = np.concatenate((TotlPoint_12to24Z_satdy_HC02,extra_X_dict["extra_HC02"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC03 = np.concatenate((TotlPoint_12to24Z_satdy_HC03,extra_X_dict["extra_HC03"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC04 = np.concatenate((TotlPoint_12to24Z_satdy_HC04,extra_X_dict["extra_HC04"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC05 = np.concatenate((TotlPoint_12to24Z_satdy_HC05,extra_X_dict["extra_HC05"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC06 = np.concatenate((TotlPoint_12to24Z_satdy_HC06,extra_X_dict["extra_HC06"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC07 = np.concatenate((TotlPoint_12to24Z_satdy_HC07,extra_X_dict["extra_HC07"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC08 = np.concatenate((TotlPoint_12to24Z_satdy_HC08,extra_X_dict["extra_HC08"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC09 = np.concatenate((TotlPoint_12to24Z_satdy_HC09,extra_X_dict["extra_HC09"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC10 = np.concatenate((TotlPoint_12to24Z_satdy_HC10,extra_X_dict["extra_HC10"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC11 = np.concatenate((TotlPoint_12to24Z_satdy_HC11,extra_X_dict["extra_HC11"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC12 = np.concatenate((TotlPoint_12to24Z_satdy_HC12,extra_X_dict["extra_HC12"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC13 = np.concatenate((TotlPoint_12to24Z_satdy_HC13,extra_X_dict["extra_HC13"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC14 = np.concatenate((TotlPoint_12to24Z_satdy_HC14,extra_X_dict["extra_HC14"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC15 = np.concatenate((TotlPoint_12to24Z_satdy_HC15,extra_X_dict["extra_HC15"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC16 = np.concatenate((TotlPoint_12to24Z_satdy_HC16,extra_X_dict["extra_HC16"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC17 = np.concatenate((TotlPoint_12to24Z_satdy_HC17,extra_X_dict["extra_HC17"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC18 = np.concatenate((TotlPoint_12to24Z_satdy_HC18,extra_X_dict["extra_HC18"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC19 = np.concatenate((TotlPoint_12to24Z_satdy_HC19,extra_X_dict["extra_HC19"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC20 = np.concatenate((TotlPoint_12to24Z_satdy_HC20,extra_X_dict["extra_HC20"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC21 = np.concatenate((TotlPoint_12to24Z_satdy_HC21,extra_X_dict["extra_HC21"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC22 = np.concatenate((TotlPoint_12to24Z_satdy_HC22,extra_X_dict["extra_HC22"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC23 = np.concatenate((TotlPoint_12to24Z_satdy_HC23,extra_X_dict["extra_HC23"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC24 = np.concatenate((TotlPoint_12to24Z_satdy_HC24,extra_X_dict["extra_HC24"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC25 = np.concatenate((TotlPoint_12to24Z_satdy_HC25,extra_X_dict["extra_HC25"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC26 = np.concatenate((TotlPoint_12to24Z_satdy_HC26,extra_X_dict["extra_HC26"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC27 = np.concatenate((TotlPoint_12to24Z_satdy_HC27,extra_X_dict["extra_HC27"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC28 = np.concatenate((TotlPoint_12to24Z_satdy_HC28,extra_X_dict["extra_HC28"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC29 = np.concatenate((TotlPoint_12to24Z_satdy_HC29,extra_X_dict["extra_HC29"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC30 = np.concatenate((TotlPoint_12to24Z_satdy_HC30,extra_X_dict["extra_HC30"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC31 = np.concatenate((TotlPoint_12to24Z_satdy_HC31,extra_X_dict["extra_HC31"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC32 = np.concatenate((TotlPoint_12to24Z_satdy_HC32,extra_X_dict["extra_HC32"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC33 = np.concatenate((TotlPoint_12to24Z_satdy_HC33,extra_X_dict["extra_HC33"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC34 = np.concatenate((TotlPoint_12to24Z_satdy_HC34,extra_X_dict["extra_HC34"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC35 = np.concatenate((TotlPoint_12to24Z_satdy_HC35,extra_X_dict["extra_HC35"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC36 = np.concatenate((TotlPoint_12to24Z_satdy_HC36,extra_X_dict["extra_HC36"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC37 = np.concatenate((TotlPoint_12to24Z_satdy_HC37,extra_X_dict["extra_HC37"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC38 = np.concatenate((TotlPoint_12to24Z_satdy_HC38,extra_X_dict["extra_HC38"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC39 = np.concatenate((TotlPoint_12to24Z_satdy_HC39,extra_X_dict["extra_HC39"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC40 = np.concatenate((TotlPoint_12to24Z_satdy_HC40,extra_X_dict["extra_HC40"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC41 = np.concatenate((TotlPoint_12to24Z_satdy_HC41,extra_X_dict["extra_HC41"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC42 = np.concatenate((TotlPoint_12to24Z_satdy_HC42,extra_X_dict["extra_HC42"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC43 = np.concatenate((TotlPoint_12to24Z_satdy_HC43,extra_X_dict["extra_HC43"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC44 = np.concatenate((TotlPoint_12to24Z_satdy_HC44,extra_X_dict["extra_HC44"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC45 = np.concatenate((TotlPoint_12to24Z_satdy_HC45,extra_X_dict["extra_HC45"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC46 = np.concatenate((TotlPoint_12to24Z_satdy_HC46,extra_X_dict["extra_HC46"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC47 = np.concatenate((TotlPoint_12to24Z_satdy_HC47,extra_X_dict["extra_HC47"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC48 = np.concatenate((TotlPoint_12to24Z_satdy_HC48,extra_X_dict["extra_HC48"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC49 = np.concatenate((TotlPoint_12to24Z_satdy_HC49,extra_X_dict["extra_HC49"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC50 = np.concatenate((TotlPoint_12to24Z_satdy_HC50,extra_X_dict["extra_HC50"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC51 = np.concatenate((TotlPoint_12to24Z_satdy_HC51,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC52 = np.concatenate((TotlPoint_12to24Z_satdy_HC52,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC53 = np.concatenate((TotlPoint_12to24Z_satdy_HC53,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC54 = np.concatenate((TotlPoint_12to24Z_satdy_HC54,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC55 = np.concatenate((TotlPoint_12to24Z_satdy_HC55,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC56 = np.concatenate((TotlPoint_12to24Z_satdy_HC56,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC57 = np.concatenate((TotlPoint_12to24Z_satdy_HC57,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC58 = np.concatenate((TotlPoint_12to24Z_satdy_HC58,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC59 = np.concatenate((TotlPoint_12to24Z_satdy_HC59,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC60 = np.concatenate((TotlPoint_12to24Z_satdy_HC60,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC61 = np.concatenate((TotlPoint_12to24Z_satdy_HC61,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC62 = np.concatenate((TotlPoint_12to24Z_satdy_HC62,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC63 = np.concatenate((TotlPoint_12to24Z_satdy_HC63,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC64 = np.concatenate((TotlPoint_12to24Z_satdy_HC64,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC65 = np.concatenate((TotlPoint_12to24Z_satdy_HC65,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC66 = np.concatenate((TotlPoint_12to24Z_satdy_HC66,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC67 = np.concatenate((TotlPoint_12to24Z_satdy_HC67,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_HC68 = np.concatenate((TotlPoint_12to24Z_satdy_HC68,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM01 = np.concatenate((TotlPoint_12to24Z_satdy_PM01,extra_X_dict["extra_PM01"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM02 = np.concatenate((TotlPoint_12to24Z_satdy_PM02,extra_X_dict["extra_PM02"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM03 = np.concatenate((TotlPoint_12to24Z_satdy_PM03,extra_X_dict["extra_PM03"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM04 = np.concatenate((TotlPoint_12to24Z_satdy_PM04,extra_X_dict["extra_PM04"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM05 = np.concatenate((TotlPoint_12to24Z_satdy_PM05,extra_X_dict["extra_PM05"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM06 = np.concatenate((TotlPoint_12to24Z_satdy_PM06,extra_X_dict["extra_PM06"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM07 = np.concatenate((TotlPoint_12to24Z_satdy_PM07,extra_X_dict["extra_PM07"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM08 = np.concatenate((TotlPoint_12to24Z_satdy_PM08,extra_X_dict["extra_PM08"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM09 = np.concatenate((TotlPoint_12to24Z_satdy_PM09,extra_X_dict["extra_PM09"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM10 = np.concatenate((TotlPoint_12to24Z_satdy_PM10,extra_X_dict["extra_PM10"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM11 = np.concatenate((TotlPoint_12to24Z_satdy_PM11,extra_X_dict["extra_PM11"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM12 = np.concatenate((TotlPoint_12to24Z_satdy_PM12,extra_X_dict["extra_PM12"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM13 = np.concatenate((TotlPoint_12to24Z_satdy_PM13,extra_X_dict["extra_PM13"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM14 = np.concatenate((TotlPoint_12to24Z_satdy_PM14,extra_X_dict["extra_PM14"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM15 = np.concatenate((TotlPoint_12to24Z_satdy_PM15,extra_X_dict["extra_PM15"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM16 = np.concatenate((TotlPoint_12to24Z_satdy_PM16,extra_X_dict["extra_PM16"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM17 = np.concatenate((TotlPoint_12to24Z_satdy_PM17,extra_X_dict["extra_PM17"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM18 = np.concatenate((TotlPoint_12to24Z_satdy_PM18,extra_X_dict["extra_PM18"]),axis=1)
TotlPoint_w_extra_12to24Z_satdy_PM19 = np.concatenate((TotlPoint_12to24Z_satdy_PM19,extra_X_dict["extra_PM19"]),axis=1)

###################################################################################################
#write total points with extra points appended
TotlPoint_w_extra_12to24Z_satdy_fn = append_dir+'/satdy/PtOG_12to24Z.nc'
TotlPoint_w_extra_12to24Z_satdy_file = Dataset(TotlPoint_w_extra_12to24Z_satdy_fn,mode='w',format='NETCDF3_64BIT')

#Creat dimensions
TotlPoint_w_extra_12to24Z_satdy_file.createDimension("ROW", nROW)
TotlPoint_w_extra_12to24Z_satdy_file.createDimension("Time", 12)
TotlPoint_w_extra_12to24Z_satdy_file.sync()

#Create variables
#float ITYPE(ROW) ;
ITYPE = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('ITYPE','f4',('ROW'),fill_value = float(0))
ITYPE[:] = TotlPoint_w_extra_12to24Z_satdy_ITYPE
varattrs=["FieldType","MemoryOrder","description","units","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['ITYPE'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['ITYPE'], varattr);
        setattr(ITYPE, varattr, varattrVal)

#float STKht(ROW) ;
STKht = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('STKht','f4',('ROW'),fill_value = float(0))
STKht[:] = TotlPoint_w_extra_12to24Z_satdy_STKht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['STKht'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['STKht'], varattr);
        setattr(STKht, varattr, varattrVal)
        
#float STKdiam(ROW) ;
STKdiam = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('STKdiam','f4',('ROW'),fill_value = float(0))
STKdiam[:] = TotlPoint_w_extra_12to24Z_satdy_STKdiam
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['STKdiam'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['STKdiam'], varattr);
        setattr(STKdiam, varattr, varattrVal)

#float STKtemp(ROW) ;
STKtemp = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('STKtemp','f4',('ROW'),fill_value = float(0))
STKtemp[:] = TotlPoint_w_extra_12to24Z_satdy_STKtemp
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['STKtemp'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['STKtemp'], varattr);
        setattr(STKtemp, varattr, varattrVal)
        
#float STKve(ROW) ;
STKve = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('STKve','f4',('ROW'),fill_value = float(0))
STKve[:] = TotlPoint_w_extra_12to24Z_satdy_STKve
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['STKve'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['STKve'], varattr);
        setattr(STKve, varattr, varattrVal)
        
#float STKflw(ROW) ;
STKflw = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('STKflw','f4',('ROW'),fill_value = float(0))
STKflw[:] = TotlPoint_w_extra_12to24Z_satdy_STKflw
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['STKflw'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['STKflw'], varattr);
        setattr(STKflw, varattr, varattrVal)
        
#float FUGht(ROW) ;
FUGht = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('FUGht','f4',('ROW'),fill_value = float(0))
FUGht[:] = TotlPoint_w_extra_12to24Z_satdy_FUGht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['FUGht'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['FUGht'], varattr);
        setattr(FUGht, varattr, varattrVal)
        
#float XLONG(ROW) ;
XLONG = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('XLONG','f4',('ROW'),fill_value = float(0))
XLONG[:] = TotlPoint_w_extra_12to24Z_satdy_XLONG
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['XLONG'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['XLONG'], varattr);
        setattr(XLONG, varattr, varattrVal)
        
#float XLAT(ROW) ;
XLAT = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('XLAT','f4',('ROW'),fill_value = float(0))
XLAT[:] = TotlPoint_w_extra_12to24Z_satdy_XLAT
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['XLAT'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['XLAT'], varattr);
        setattr(XLAT, varattr, varattrVal)

#float CO2(Time, ROW) ;
CO2 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('CO2','f4',('Time','ROW'))
CO2[:] = TotlPoint_w_extra_12to24Z_satdy_CO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['CO2'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['CO2'], varattr);
        setattr(CO2, varattr, varattrVal)
        
#float CO(Time, ROW) ;
CO = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('CO','f4',('Time','ROW'))
CO[:] = TotlPoint_w_extra_12to24Z_satdy_CO
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['CO'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['CO'], varattr);
        setattr(CO, varattr, varattrVal)
        
#float NH3(Time, ROW) ;
NH3 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('NH3','f4',('Time','ROW'))
NH3[:] = TotlPoint_w_extra_12to24Z_satdy_NH3
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['NH3'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['NH3'], varattr);
        setattr(NH3, varattr, varattrVal)
        
#float NOX(Time, ROW) ;
NOX = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('NOX','f4',('Time','ROW'))
NOX[:] = TotlPoint_w_extra_12to24Z_satdy_NOX
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['NOX'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['NOX'], varattr);
        setattr(NOX, varattr, varattrVal)
        
#float PM10-PRI(Time, ROW) ;
PM10_PRI = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM10-PRI','f4',('Time','ROW'))
PM10_PRI[:] = TotlPoint_w_extra_12to24Z_satdy_PM10_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM10-PRI'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM10-PRI'], varattr);
        setattr(PM10_PRI, varattr, varattrVal)
        
#float PM25-PRI(Time, ROW) ;
PM25_PRI = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM25-PRI','f4',('Time','ROW'))
PM25_PRI[:] = TotlPoint_w_extra_12to24Z_satdy_PM25_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM25-PRI'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM25-PRI'], varattr);
        setattr(PM25_PRI, varattr, varattrVal)
        
#float SO2(Time, ROW) ;
SO2 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('SO2','f4',('Time','ROW'))
SO2[:] = TotlPoint_w_extra_12to24Z_satdy_SO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['SO2'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['SO2'], varattr);
        setattr(SO2, varattr, varattrVal)
        
#float VOC(Time, ROW) ;
VOC = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('VOC','f4',('Time','ROW'))
VOC[:] = TotlPoint_w_extra_12to24Z_satdy_VOC
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['VOC'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['VOC'], varattr);
        setattr(VOC, varattr, varattrVal)
        
#float HC01(Time, ROW) ;
HC01 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC01','f4',('Time','ROW'))
HC01[:] = TotlPoint_w_extra_12to24Z_satdy_HC01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC01'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC01'], varattr);
        setattr(HC01, varattr, varattrVal)
        
#float HC02(Time, ROW) ;
HC02 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC02','f4',('Time','ROW'))
HC02[:] = TotlPoint_w_extra_12to24Z_satdy_HC02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC02'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC02'], varattr);
        setattr(HC02, varattr, varattrVal)
        
#float HC03(Time, ROW) ;
HC03 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC03','f4',('Time','ROW'))
HC03[:] = TotlPoint_w_extra_12to24Z_satdy_HC03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC03'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC03'], varattr);
        setattr(HC03, varattr, varattrVal)
        
#float HC04(Time, ROW) ;
HC04 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC04','f4',('Time','ROW'))
HC04[:] = TotlPoint_w_extra_12to24Z_satdy_HC04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC04'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC04'], varattr);
        setattr(HC04, varattr, varattrVal)
        
#float HC05(Time, ROW) ;
HC05 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC05','f4',('Time','ROW'))
HC05[:] = TotlPoint_w_extra_12to24Z_satdy_HC05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC05'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC05'], varattr);
        setattr(HC05, varattr, varattrVal)
        
#float HC06(Time, ROW) ;
HC06 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC06','f4',('Time','ROW'))
HC06[:] = TotlPoint_w_extra_12to24Z_satdy_HC06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC06'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC06'], varattr);
        setattr(HC06, varattr, varattrVal)
        
#float HC07(Time, ROW) ;
HC07 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC07','f4',('Time','ROW'))
HC07[:] = TotlPoint_w_extra_12to24Z_satdy_HC07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC07'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC07'], varattr);
        setattr(HC07, varattr, varattrVal)
        
#float HC08(Time, ROW) ;
HC08 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC08','f4',('Time','ROW'))
HC08[:] = TotlPoint_w_extra_12to24Z_satdy_HC08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC08'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC08'], varattr);
        setattr(HC08, varattr, varattrVal)
        
#float HC09(Time, ROW) ;
HC09 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC09','f4',('Time','ROW'))
HC09[:] = TotlPoint_w_extra_12to24Z_satdy_HC09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC09'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC09'], varattr);
        setattr(HC09, varattr, varattrVal)
        
#float HC10(Time, ROW) ;
HC10 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC10','f4',('Time','ROW'))
HC10[:] = TotlPoint_w_extra_12to24Z_satdy_HC10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC10'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC10'], varattr);
        setattr(HC10, varattr, varattrVal)

#float HC11(Time, ROW) ;
HC11 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC11','f4',('Time','ROW'))
HC11[:] = TotlPoint_w_extra_12to24Z_satdy_HC11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC11'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC11'], varattr);
        setattr(HC11, varattr, varattrVal)
        
#float HC12(Time, ROW) ;
HC12 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC12','f4',('Time','ROW'))
HC12[:] = TotlPoint_w_extra_12to24Z_satdy_HC12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC12'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC12'], varattr);
        setattr(HC12, varattr, varattrVal)
        
#float HC13(Time, ROW) ;
HC13 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC13','f4',('Time','ROW'))
HC13[:] = TotlPoint_w_extra_12to24Z_satdy_HC13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC13'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC13'], varattr);
        setattr(HC13, varattr, varattrVal)
        
#float HC14(Time, ROW) ;
HC14 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC14','f4',('Time','ROW'))
HC14[:] = TotlPoint_w_extra_12to24Z_satdy_HC14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC14'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC14'], varattr);
        setattr(HC14, varattr, varattrVal)
        
#float HC15(Time, ROW) ;
HC15 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC15','f4',('Time','ROW'))
HC15[:] = TotlPoint_w_extra_12to24Z_satdy_HC15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC15'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC15'], varattr);
        setattr(HC15, varattr, varattrVal)
        
#float HC16(Time, ROW) ;
HC16 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC16','f4',('Time','ROW'))
HC16[:] = TotlPoint_w_extra_12to24Z_satdy_HC16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC16'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC16'], varattr);
        setattr(HC16, varattr, varattrVal)
        
#float HC17(Time, ROW) ;
HC17 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC17','f4',('Time','ROW'))
HC17[:] = TotlPoint_w_extra_12to24Z_satdy_HC17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC17'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC17'], varattr);
        setattr(HC17, varattr, varattrVal)
        
#float HC18(Time, ROW) ;
HC18 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC18','f4',('Time','ROW'))
HC18[:] = TotlPoint_w_extra_12to24Z_satdy_HC18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC18'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC18'], varattr);
        setattr(HC18, varattr, varattrVal)
        
#float HC19(Time, ROW) ;
HC19 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC19','f4',('Time','ROW'))
HC19[:] = TotlPoint_w_extra_12to24Z_satdy_HC19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC19'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC19'], varattr);
        setattr(HC19, varattr, varattrVal)

#float HC20(Time, ROW) ;
HC20 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC20','f4',('Time','ROW'))
HC20[:] = TotlPoint_w_extra_12to24Z_satdy_HC20
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC20'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC20'], varattr);
        setattr(HC20, varattr, varattrVal)

#float HC21(Time, ROW) ;
HC21 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC21','f4',('Time','ROW'))
HC21[:] = TotlPoint_w_extra_12to24Z_satdy_HC21
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC21'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC21'], varattr);
        setattr(HC21, varattr, varattrVal)
        
#float HC22(Time, ROW) ;
HC22 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC22','f4',('Time','ROW'))
HC22[:] = TotlPoint_w_extra_12to24Z_satdy_HC22
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC22'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC22'], varattr);
        setattr(HC22, varattr, varattrVal)
        
#float HC23(Time, ROW) ;
HC23 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC23','f4',('Time','ROW'))
HC23[:] = TotlPoint_w_extra_12to24Z_satdy_HC23
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC23'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC23'], varattr);
        setattr(HC23, varattr, varattrVal)
        
#float HC24(Time, ROW) ;
HC24 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC24','f4',('Time','ROW'))
HC24[:] = TotlPoint_w_extra_12to24Z_satdy_HC24
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC24'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC24'], varattr);
        setattr(HC24, varattr, varattrVal)
        
#float HC25(Time, ROW) ;
HC25 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC25','f4',('Time','ROW'))
HC25[:] = TotlPoint_w_extra_12to24Z_satdy_HC25
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC25'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC25'], varattr);
        setattr(HC25, varattr, varattrVal)
        
#float HC26(Time, ROW) ;
HC26 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC26','f4',('Time','ROW'))
HC26[:] = TotlPoint_w_extra_12to24Z_satdy_HC26
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC26'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC26'], varattr);
        setattr(HC26, varattr, varattrVal)
        
#float HC27(Time, ROW) ;
HC27 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC27','f4',('Time','ROW'))
HC27[:] = TotlPoint_w_extra_12to24Z_satdy_HC27
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC27'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC27'], varattr);
        setattr(HC27, varattr, varattrVal)
        
#float HC28(Time, ROW) ;
HC28 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC28','f4',('Time','ROW'))
HC28[:] = TotlPoint_w_extra_12to24Z_satdy_HC28
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC28'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC28'], varattr);
        setattr(HC28, varattr, varattrVal)
        
#float HC29(Time, ROW) ;
HC29 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC29','f4',('Time','ROW'))
HC29[:] = TotlPoint_w_extra_12to24Z_satdy_HC29
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC29'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC29'], varattr);
        setattr(HC29, varattr, varattrVal)

#float HC30(Time, ROW) ;
HC30 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC30','f4',('Time','ROW'))
HC30[:] = TotlPoint_w_extra_12to24Z_satdy_HC30
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC30'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC30'], varattr);
        setattr(HC30, varattr, varattrVal)

#float HC31(Time, ROW) ;
HC31 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC31','f4',('Time','ROW'))
HC31[:] = TotlPoint_w_extra_12to24Z_satdy_HC31
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC31'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC31'], varattr);
        setattr(HC31, varattr, varattrVal)
        
#float HC32(Time, ROW) ;
HC32 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC32','f4',('Time','ROW'))
HC32[:] = TotlPoint_w_extra_12to24Z_satdy_HC32
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC32'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC32'], varattr);
        setattr(HC32, varattr, varattrVal)
        
#float HC33(Time, ROW) ;
HC33 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC33','f4',('Time','ROW'))
HC33[:] = TotlPoint_w_extra_12to24Z_satdy_HC33
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC33'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC33'], varattr);
        setattr(HC33, varattr, varattrVal)
        
#float HC34(Time, ROW) ;
HC34 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC34','f4',('Time','ROW'))
HC34[:] = TotlPoint_w_extra_12to24Z_satdy_HC34
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC34'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC34'], varattr);
        setattr(HC34, varattr, varattrVal)
        
#float HC35(Time, ROW) ;
HC35 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC35','f4',('Time','ROW'))
HC35[:] = TotlPoint_w_extra_12to24Z_satdy_HC35
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC35'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC35'], varattr);
        setattr(HC35, varattr, varattrVal)
        
#float HC36(Time, ROW) ;
HC36 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC36','f4',('Time','ROW'))
HC36[:] = TotlPoint_w_extra_12to24Z_satdy_HC36
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC36'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC36'], varattr);
        setattr(HC36, varattr, varattrVal)
        
#float HC37(Time, ROW) ;
HC37 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC37','f4',('Time','ROW'))
HC37[:] = TotlPoint_w_extra_12to24Z_satdy_HC37
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC37'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC37'], varattr);
        setattr(HC37, varattr, varattrVal)
        
#float HC38(Time, ROW) ;
HC38 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC38','f4',('Time','ROW'))
HC38[:] = TotlPoint_w_extra_12to24Z_satdy_HC38
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC38'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC38'], varattr);
        setattr(HC38, varattr, varattrVal)
        
#float HC39(Time, ROW) ;
HC39 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC39','f4',('Time','ROW'))
HC39[:] = TotlPoint_w_extra_12to24Z_satdy_HC39
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC39'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC39'], varattr);
        setattr(HC39, varattr, varattrVal)
        
#float HC40(Time, ROW) ;
HC40 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC40','f4',('Time','ROW'))
HC40[:] = TotlPoint_w_extra_12to24Z_satdy_HC40
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC40'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC40'], varattr);
        setattr(HC40, varattr, varattrVal)

#float HC41(Time, ROW) ;
HC41 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC41','f4',('Time','ROW'))
HC41[:] = TotlPoint_w_extra_12to24Z_satdy_HC41
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC41'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC41'], varattr);
        setattr(HC41, varattr, varattrVal)
        
#float HC42(Time, ROW) ;
HC42 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC42','f4',('Time','ROW'))
HC42[:] = TotlPoint_w_extra_12to24Z_satdy_HC42
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC42'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC42'], varattr);
        setattr(HC42, varattr, varattrVal)
        
#float HC43(Time, ROW) ;
HC43 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC43','f4',('Time','ROW'))
HC43[:] = TotlPoint_w_extra_12to24Z_satdy_HC43
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC43'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC43'], varattr);
        setattr(HC43, varattr, varattrVal)
        
#float HC44(Time, ROW) ;
HC44 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC44','f4',('Time','ROW'))
HC44[:] = TotlPoint_w_extra_12to24Z_satdy_HC44
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC44'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC44'], varattr);
        setattr(HC44, varattr, varattrVal)
        
#float HC45(Time, ROW) ;
HC45 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC45','f4',('Time','ROW'))
HC45[:] = TotlPoint_w_extra_12to24Z_satdy_HC45
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC45'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC45'], varattr);
        setattr(HC45, varattr, varattrVal)
        
#float HC46(Time, ROW) ;
HC46 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC46','f4',('Time','ROW'))
HC46[:] = TotlPoint_w_extra_12to24Z_satdy_HC46
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC46'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC46'], varattr);
        setattr(HC46, varattr, varattrVal)
        
#float HC47(Time, ROW) ;
HC47 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC47','f4',('Time','ROW'))
HC47[:] = TotlPoint_w_extra_12to24Z_satdy_HC47
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC47'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC47'], varattr);
        setattr(HC47, varattr, varattrVal)
        
#float HC48(Time, ROW) ;
HC48 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC48','f4',('Time','ROW'))
HC48[:] = TotlPoint_w_extra_12to24Z_satdy_HC48
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC48'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC48'], varattr);
        setattr(HC48, varattr, varattrVal)
        
#float HC49(Time, ROW) ;
HC49 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC49','f4',('Time','ROW'))
HC49[:] = TotlPoint_w_extra_12to24Z_satdy_HC49
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC49'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC49'], varattr);
        setattr(HC49, varattr, varattrVal)
        
#float HC50(Time, ROW) ;
HC50 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC50','f4',('Time','ROW'))
HC50[:] = TotlPoint_w_extra_12to24Z_satdy_HC50
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC50'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC50'], varattr);
        setattr(HC50, varattr, varattrVal)

#float HC51(Time, ROW) ;
HC51 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC51','f4',('Time','ROW'))
HC51[:] = TotlPoint_w_extra_12to24Z_satdy_HC51
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC51'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC51'], varattr);
        setattr(HC51, varattr, varattrVal)
        
#float HC52(Time, ROW) ;
HC52 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC52','f4',('Time','ROW'))
HC52[:] = TotlPoint_w_extra_12to24Z_satdy_HC52
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC52'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC52'], varattr);
        setattr(HC52, varattr, varattrVal)
        
#float HC53(Time, ROW) ;
HC53 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC53','f4',('Time','ROW'))
HC53[:] = TotlPoint_w_extra_12to24Z_satdy_HC53
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC53'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC53'], varattr);
        setattr(HC53, varattr, varattrVal)
        
#float HC54(Time, ROW) ;
HC54 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC54','f4',('Time','ROW'))
HC54[:] = TotlPoint_w_extra_12to24Z_satdy_HC54
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC54'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC54'], varattr);
        setattr(HC54, varattr, varattrVal)
        
#float HC55(Time, ROW) ;
HC55 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC55','f4',('Time','ROW'))
HC55[:] = TotlPoint_w_extra_12to24Z_satdy_HC55
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC55'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC55'], varattr);
        setattr(HC55, varattr, varattrVal)
        
#float HC56(Time, ROW) ;
HC56 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC56','f4',('Time','ROW'))
HC56[:] = TotlPoint_w_extra_12to24Z_satdy_HC56
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC56'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC56'], varattr);
        setattr(HC56, varattr, varattrVal)
        
#float HC57(Time, ROW) ;
HC57 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC57','f4',('Time','ROW'))
HC57[:] = TotlPoint_w_extra_12to24Z_satdy_HC57
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC57'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC57'], varattr);
        setattr(HC57, varattr, varattrVal)
        
#float HC58(Time, ROW) ;
HC58 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC58','f4',('Time','ROW'))
HC58[:] = TotlPoint_w_extra_12to24Z_satdy_HC58
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC58'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC58'], varattr);
        setattr(HC58, varattr, varattrVal)
        
#float HC59(Time, ROW) ;
HC59 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC59','f4',('Time','ROW'))
HC59[:] = TotlPoint_w_extra_12to24Z_satdy_HC59
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC59'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC59'], varattr);
        setattr(HC59, varattr, varattrVal)

#float HC60(Time, ROW) ;
HC60 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC60','f4',('Time','ROW'))
HC60[:] = TotlPoint_w_extra_12to24Z_satdy_HC60
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC60'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC60'], varattr);
        setattr(HC60, varattr, varattrVal)

#float HC61(Time, ROW) ;
HC61 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC61','f4',('Time','ROW'))
HC61[:] = TotlPoint_w_extra_12to24Z_satdy_HC61
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC61'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC61'], varattr);
        setattr(HC61, varattr, varattrVal)
        
#float HC62(Time, ROW) ;
HC62 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC62','f4',('Time','ROW'))
HC62[:] = TotlPoint_w_extra_12to24Z_satdy_HC62
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC62'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC62'], varattr);
        setattr(HC62, varattr, varattrVal)
        
#float HC63(Time, ROW) ;
HC63 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC63','f4',('Time','ROW'))
HC63[:] = TotlPoint_w_extra_12to24Z_satdy_HC63
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC63'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC63'], varattr);
        setattr(HC63, varattr, varattrVal)
        
#float HC64(Time, ROW) ;
HC64 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC64','f4',('Time','ROW'))
HC64[:] = TotlPoint_w_extra_12to24Z_satdy_HC64
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC64'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC64'], varattr);
        setattr(HC64, varattr, varattrVal)
        
#float HC65(Time, ROW) ;
HC65 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC65','f4',('Time','ROW'))
HC65[:] = TotlPoint_w_extra_12to24Z_satdy_HC65
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC65'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC65'], varattr);
        setattr(HC65, varattr, varattrVal)
        
#float HC66(Time, ROW) ;
HC66 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC66','f4',('Time','ROW'))
HC66[:] = TotlPoint_w_extra_12to24Z_satdy_HC66
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC66'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC66'], varattr);
        setattr(HC66, varattr, varattrVal)
        
#float HC67(Time, ROW) ;
HC67 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC67','f4',('Time','ROW'))
HC67[:] = TotlPoint_w_extra_12to24Z_satdy_HC67
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC67'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC67'], varattr);
        setattr(HC67, varattr, varattrVal)
        
#float HC68(Time, ROW) ;
HC68 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('HC68','f4',('Time','ROW'))
HC68[:] = TotlPoint_w_extra_12to24Z_satdy_HC68
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['HC68'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['HC68'], varattr);
        setattr(HC68, varattr, varattrVal)

#float PM01(Time, ROW) ;
PM01 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM01','f4',('Time','ROW'))
PM01[:] = TotlPoint_w_extra_12to24Z_satdy_PM01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM01'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM01'], varattr);
        setattr(PM01, varattr, varattrVal)

#float PM02(Time, ROW) ;
PM02 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM02','f4',('Time','ROW'))
PM02[:] = TotlPoint_w_extra_12to24Z_satdy_PM02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM02'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM02'], varattr);
        setattr(PM02, varattr, varattrVal)
        
#float PM03(Time, ROW) ;
PM03 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM03','f4',('Time','ROW'))
PM03[:] = TotlPoint_w_extra_12to24Z_satdy_PM03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM03'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM03'], varattr);
        setattr(PM03, varattr, varattrVal)

#float PM04(Time, ROW) ;
PM04 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM04','f4',('Time','ROW'))
PM04[:] = TotlPoint_w_extra_12to24Z_satdy_PM04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM04'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM04'], varattr);
        setattr(PM04, varattr, varattrVal)
        
#float PM05(Time, ROW) ;
PM05 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM05','f4',('Time','ROW'))
PM05[:] = TotlPoint_w_extra_12to24Z_satdy_PM05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM05'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM05'], varattr);
        setattr(PM05, varattr, varattrVal)

#float PM06(Time, ROW) ;
PM06 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM06','f4',('Time','ROW'))
PM06[:] = TotlPoint_w_extra_12to24Z_satdy_PM06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM06'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM06'], varattr);
        setattr(PM06, varattr, varattrVal)
        
#float PM07(Time, ROW) ;
PM07 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM07','f4',('Time','ROW'))
PM07[:] = TotlPoint_w_extra_12to24Z_satdy_PM07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM07'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM07'], varattr);
        setattr(PM07, varattr, varattrVal)
        
#float PM08(Time, ROW) ;
PM08 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM08','f4',('Time','ROW'))
PM08[:] = TotlPoint_w_extra_12to24Z_satdy_PM08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM08'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM08'], varattr);
        setattr(PM08, varattr, varattrVal)
        
#float PM09(Time, ROW) ;
PM09 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM09','f4',('Time','ROW'))
PM09[:] = TotlPoint_w_extra_12to24Z_satdy_PM09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM09'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM09'], varattr);
        setattr(PM09, varattr, varattrVal)
        
#float PM10(Time, ROW) ;
PM10 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM10','f4',('Time','ROW'))
PM10[:] = TotlPoint_w_extra_12to24Z_satdy_PM10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM10'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM10'], varattr);
        setattr(PM10, varattr, varattrVal)
        
#float PM11(Time, ROW) ;
PM11 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM11','f4',('Time','ROW'))
PM11[:] = TotlPoint_w_extra_12to24Z_satdy_PM11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM11'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM11'], varattr);
        setattr(PM11, varattr, varattrVal)

#float PM12(Time, ROW) ;
PM12 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM12','f4',('Time','ROW'))
PM12[:] = TotlPoint_w_extra_12to24Z_satdy_PM12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM12'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM12'], varattr);
        setattr(PM12, varattr, varattrVal)
        
#float PM13(Time, ROW) ;
PM13 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM13','f4',('Time','ROW'))
PM13[:] = TotlPoint_w_extra_12to24Z_satdy_PM13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM13'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM13'], varattr);
        setattr(PM13, varattr, varattrVal)

#float PM14(Time, ROW) ;
PM14 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM14','f4',('Time','ROW'))
PM14[:] = TotlPoint_w_extra_12to24Z_satdy_PM14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM14'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM14'], varattr);
        setattr(PM14, varattr, varattrVal)
        
#float PM15(Time, ROW) ;
PM15 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM15','f4',('Time','ROW'))
PM15[:] = TotlPoint_w_extra_12to24Z_satdy_PM15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM15'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM15'], varattr);
        setattr(PM15, varattr, varattrVal)

#float PM16(Time, ROW) ;
PM16 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM16','f4',('Time','ROW'))
PM16[:] = TotlPoint_w_extra_12to24Z_satdy_PM16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM16'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM16'], varattr);
        setattr(PM16, varattr, varattrVal)
        
#float PM17(Time, ROW) ;
PM17 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM17','f4',('Time','ROW'))
PM17[:] = TotlPoint_w_extra_12to24Z_satdy_PM17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM17'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM17'], varattr);
        setattr(PM17, varattr, varattrVal)
        
#float PM18(Time, ROW) ;
PM18 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM18','f4',('Time','ROW'))
PM18[:] = TotlPoint_w_extra_12to24Z_satdy_PM18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM18'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM18'], varattr);
        setattr(PM18, varattr, varattrVal)
        
#float PM19(Time, ROW) ;
PM19 = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('PM19','f4',('Time','ROW'))
PM19[:] = TotlPoint_w_extra_12to24Z_satdy_PM19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_satdy_file.variables['PM19'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file.variables['PM19'], varattr);
        setattr(PM19, varattr, varattrVal)

#char Times(Time) ;
Times = TotlPoint_w_extra_12to24Z_satdy_file.createVariable('Times','S1',('Time'))
Times[:] = TotlPoint_12to24Z_satdy_Times

#copy global attributes from TotlPoint_12to24Z_satdy_file
for varattr in TotlPoint_12to24Z_satdy_file.ncattrs():
    if hasattr(TotlPoint_12to24Z_satdy_file, varattr):
        varattrVal = getattr(TotlPoint_12to24Z_satdy_file, varattr);
        setattr(TotlPoint_w_extra_12to24Z_satdy_file, varattr, varattrVal)

TotlPoint_w_extra_12to24Z_satdy_file.close()


# In[20]:


#append extra points to original point file that is input to wrfchemi assembly program

###################################################################################################
#sundy, 00to12Z

###################################################################################################
#read original variables
TotlPoint_00to12Z_sundy_fn = base_dir+'/sundy/PtOG_00to12Z.nc'
TotlPoint_00to12Z_sundy_file = Dataset(TotlPoint_00to12Z_sundy_fn,mode='r',open=True)
TotlPoint_00to12Z_sundy_ITYPE = TotlPoint_00to12Z_sundy_file.variables['ITYPE'][:]
TotlPoint_00to12Z_sundy_STKht = TotlPoint_00to12Z_sundy_file.variables['STKht'][:]
TotlPoint_00to12Z_sundy_STKdiam = TotlPoint_00to12Z_sundy_file.variables['STKdiam'][:]
TotlPoint_00to12Z_sundy_STKtemp = TotlPoint_00to12Z_sundy_file.variables['STKtemp'][:]
TotlPoint_00to12Z_sundy_STKve = TotlPoint_00to12Z_sundy_file.variables['STKve'][:]
TotlPoint_00to12Z_sundy_STKflw = TotlPoint_00to12Z_sundy_file.variables['STKflw'][:]
TotlPoint_00to12Z_sundy_FUGht = TotlPoint_00to12Z_sundy_file.variables['FUGht'][:]
TotlPoint_00to12Z_sundy_XLONG = TotlPoint_00to12Z_sundy_file.variables['XLONG'][:]
TotlPoint_00to12Z_sundy_XLAT = TotlPoint_00to12Z_sundy_file.variables['XLAT'][:]
TotlPoint_00to12Z_sundy_CO2 = TotlPoint_00to12Z_sundy_file.variables['CO2'][:][:] 
TotlPoint_00to12Z_sundy_CO = TotlPoint_00to12Z_sundy_file.variables['CO'][:][:] 
TotlPoint_00to12Z_sundy_NH3 = TotlPoint_00to12Z_sundy_file.variables['NH3'][:][:] 
TotlPoint_00to12Z_sundy_NOX = TotlPoint_00to12Z_sundy_file.variables['NOX'][:][:] 
TotlPoint_00to12Z_sundy_PM10_PRI = TotlPoint_00to12Z_sundy_file.variables['PM10-PRI'][:][:] 
TotlPoint_00to12Z_sundy_PM25_PRI = TotlPoint_00to12Z_sundy_file.variables['PM25-PRI'][:][:] 
TotlPoint_00to12Z_sundy_SO2 = TotlPoint_00to12Z_sundy_file.variables['SO2'][:][:] 
TotlPoint_00to12Z_sundy_VOC = TotlPoint_00to12Z_sundy_file.variables['VOC'][:][:] 
TotlPoint_00to12Z_sundy_HC01 = TotlPoint_00to12Z_sundy_file.variables['HC01'][:][:] 
TotlPoint_00to12Z_sundy_HC02 = TotlPoint_00to12Z_sundy_file.variables['HC02'][:][:] 
TotlPoint_00to12Z_sundy_HC03 = TotlPoint_00to12Z_sundy_file.variables['HC03'][:][:] 
TotlPoint_00to12Z_sundy_HC04 = TotlPoint_00to12Z_sundy_file.variables['HC04'][:][:] 
TotlPoint_00to12Z_sundy_HC05 = TotlPoint_00to12Z_sundy_file.variables['HC05'][:][:] 
TotlPoint_00to12Z_sundy_HC06 = TotlPoint_00to12Z_sundy_file.variables['HC06'][:][:] 
TotlPoint_00to12Z_sundy_HC07 = TotlPoint_00to12Z_sundy_file.variables['HC07'][:][:] 
TotlPoint_00to12Z_sundy_HC08 = TotlPoint_00to12Z_sundy_file.variables['HC08'][:][:] 
TotlPoint_00to12Z_sundy_HC09 = TotlPoint_00to12Z_sundy_file.variables['HC09'][:][:] 
TotlPoint_00to12Z_sundy_HC10 = TotlPoint_00to12Z_sundy_file.variables['HC10'][:][:] 
TotlPoint_00to12Z_sundy_HC11 = TotlPoint_00to12Z_sundy_file.variables['HC11'][:][:] 
TotlPoint_00to12Z_sundy_HC12 = TotlPoint_00to12Z_sundy_file.variables['HC12'][:][:] 
TotlPoint_00to12Z_sundy_HC13 = TotlPoint_00to12Z_sundy_file.variables['HC13'][:][:] 
TotlPoint_00to12Z_sundy_HC14 = TotlPoint_00to12Z_sundy_file.variables['HC14'][:][:] 
TotlPoint_00to12Z_sundy_HC15 = TotlPoint_00to12Z_sundy_file.variables['HC15'][:][:] 
TotlPoint_00to12Z_sundy_HC16 = TotlPoint_00to12Z_sundy_file.variables['HC16'][:][:] 
TotlPoint_00to12Z_sundy_HC17 = TotlPoint_00to12Z_sundy_file.variables['HC17'][:][:] 
TotlPoint_00to12Z_sundy_HC18 = TotlPoint_00to12Z_sundy_file.variables['HC18'][:][:] 
TotlPoint_00to12Z_sundy_HC19 = TotlPoint_00to12Z_sundy_file.variables['HC19'][:][:] 
TotlPoint_00to12Z_sundy_HC20 = TotlPoint_00to12Z_sundy_file.variables['HC20'][:][:] 
TotlPoint_00to12Z_sundy_HC21 = TotlPoint_00to12Z_sundy_file.variables['HC21'][:][:] 
TotlPoint_00to12Z_sundy_HC22 = TotlPoint_00to12Z_sundy_file.variables['HC22'][:][:] 
TotlPoint_00to12Z_sundy_HC23 = TotlPoint_00to12Z_sundy_file.variables['HC23'][:][:] 
TotlPoint_00to12Z_sundy_HC24 = TotlPoint_00to12Z_sundy_file.variables['HC24'][:][:] 
TotlPoint_00to12Z_sundy_HC25 = TotlPoint_00to12Z_sundy_file.variables['HC25'][:][:] 
TotlPoint_00to12Z_sundy_HC26 = TotlPoint_00to12Z_sundy_file.variables['HC26'][:][:] 
TotlPoint_00to12Z_sundy_HC27 = TotlPoint_00to12Z_sundy_file.variables['HC27'][:][:] 
TotlPoint_00to12Z_sundy_HC28 = TotlPoint_00to12Z_sundy_file.variables['HC28'][:][:] 
TotlPoint_00to12Z_sundy_HC29 = TotlPoint_00to12Z_sundy_file.variables['HC29'][:][:] 
TotlPoint_00to12Z_sundy_HC30 = TotlPoint_00to12Z_sundy_file.variables['HC30'][:][:] 
TotlPoint_00to12Z_sundy_HC31 = TotlPoint_00to12Z_sundy_file.variables['HC31'][:][:] 
TotlPoint_00to12Z_sundy_HC32 = TotlPoint_00to12Z_sundy_file.variables['HC32'][:][:] 
TotlPoint_00to12Z_sundy_HC33 = TotlPoint_00to12Z_sundy_file.variables['HC33'][:][:] 
TotlPoint_00to12Z_sundy_HC34 = TotlPoint_00to12Z_sundy_file.variables['HC34'][:][:] 
TotlPoint_00to12Z_sundy_HC35 = TotlPoint_00to12Z_sundy_file.variables['HC35'][:][:] 
TotlPoint_00to12Z_sundy_HC36 = TotlPoint_00to12Z_sundy_file.variables['HC36'][:][:] 
TotlPoint_00to12Z_sundy_HC37 = TotlPoint_00to12Z_sundy_file.variables['HC37'][:][:] 
TotlPoint_00to12Z_sundy_HC38 = TotlPoint_00to12Z_sundy_file.variables['HC38'][:][:] 
TotlPoint_00to12Z_sundy_HC39 = TotlPoint_00to12Z_sundy_file.variables['HC39'][:][:] 
TotlPoint_00to12Z_sundy_HC40 = TotlPoint_00to12Z_sundy_file.variables['HC40'][:][:] 
TotlPoint_00to12Z_sundy_HC41 = TotlPoint_00to12Z_sundy_file.variables['HC41'][:][:] 
TotlPoint_00to12Z_sundy_HC42 = TotlPoint_00to12Z_sundy_file.variables['HC42'][:][:] 
TotlPoint_00to12Z_sundy_HC43 = TotlPoint_00to12Z_sundy_file.variables['HC43'][:][:] 
TotlPoint_00to12Z_sundy_HC44 = TotlPoint_00to12Z_sundy_file.variables['HC44'][:][:] 
TotlPoint_00to12Z_sundy_HC45 = TotlPoint_00to12Z_sundy_file.variables['HC45'][:][:] 
TotlPoint_00to12Z_sundy_HC46 = TotlPoint_00to12Z_sundy_file.variables['HC46'][:][:] 
TotlPoint_00to12Z_sundy_HC47 = TotlPoint_00to12Z_sundy_file.variables['HC47'][:][:] 
TotlPoint_00to12Z_sundy_HC48 = TotlPoint_00to12Z_sundy_file.variables['HC48'][:][:] 
TotlPoint_00to12Z_sundy_HC49 = TotlPoint_00to12Z_sundy_file.variables['HC49'][:][:] 
TotlPoint_00to12Z_sundy_HC50 = TotlPoint_00to12Z_sundy_file.variables['HC50'][:][:] 
TotlPoint_00to12Z_sundy_HC51 = TotlPoint_00to12Z_sundy_file.variables['HC51'][:][:] 
TotlPoint_00to12Z_sundy_HC52 = TotlPoint_00to12Z_sundy_file.variables['HC52'][:][:] 
TotlPoint_00to12Z_sundy_HC53 = TotlPoint_00to12Z_sundy_file.variables['HC53'][:][:] 
TotlPoint_00to12Z_sundy_HC54 = TotlPoint_00to12Z_sundy_file.variables['HC54'][:][:] 
TotlPoint_00to12Z_sundy_HC55 = TotlPoint_00to12Z_sundy_file.variables['HC55'][:][:] 
TotlPoint_00to12Z_sundy_HC56 = TotlPoint_00to12Z_sundy_file.variables['HC56'][:][:] 
TotlPoint_00to12Z_sundy_HC57 = TotlPoint_00to12Z_sundy_file.variables['HC57'][:][:] 
TotlPoint_00to12Z_sundy_HC58 = TotlPoint_00to12Z_sundy_file.variables['HC58'][:][:] 
TotlPoint_00to12Z_sundy_HC59 = TotlPoint_00to12Z_sundy_file.variables['HC59'][:][:] 
TotlPoint_00to12Z_sundy_HC60 = TotlPoint_00to12Z_sundy_file.variables['HC60'][:][:] 
TotlPoint_00to12Z_sundy_HC61 = TotlPoint_00to12Z_sundy_file.variables['HC61'][:][:] 
TotlPoint_00to12Z_sundy_HC62 = TotlPoint_00to12Z_sundy_file.variables['HC62'][:][:] 
TotlPoint_00to12Z_sundy_HC63 = TotlPoint_00to12Z_sundy_file.variables['HC63'][:][:] 
TotlPoint_00to12Z_sundy_HC64 = TotlPoint_00to12Z_sundy_file.variables['HC64'][:][:] 
TotlPoint_00to12Z_sundy_HC65 = TotlPoint_00to12Z_sundy_file.variables['HC65'][:][:] 
TotlPoint_00to12Z_sundy_HC66 = TotlPoint_00to12Z_sundy_file.variables['HC66'][:][:] 
TotlPoint_00to12Z_sundy_HC67 = TotlPoint_00to12Z_sundy_file.variables['HC67'][:][:] 
TotlPoint_00to12Z_sundy_HC68 = TotlPoint_00to12Z_sundy_file.variables['HC68'][:][:] 
TotlPoint_00to12Z_sundy_PM01 = TotlPoint_00to12Z_sundy_file.variables['PM01'][:][:] 
TotlPoint_00to12Z_sundy_PM02 = TotlPoint_00to12Z_sundy_file.variables['PM02'][:][:] 
TotlPoint_00to12Z_sundy_PM03 = TotlPoint_00to12Z_sundy_file.variables['PM03'][:][:] 
TotlPoint_00to12Z_sundy_PM04 = TotlPoint_00to12Z_sundy_file.variables['PM04'][:][:] 
TotlPoint_00to12Z_sundy_PM05 = TotlPoint_00to12Z_sundy_file.variables['PM05'][:][:] 
TotlPoint_00to12Z_sundy_PM06 = TotlPoint_00to12Z_sundy_file.variables['PM06'][:][:] 
TotlPoint_00to12Z_sundy_PM07 = TotlPoint_00to12Z_sundy_file.variables['PM07'][:][:] 
TotlPoint_00to12Z_sundy_PM08 = TotlPoint_00to12Z_sundy_file.variables['PM08'][:][:] 
TotlPoint_00to12Z_sundy_PM09 = TotlPoint_00to12Z_sundy_file.variables['PM09'][:][:] 
TotlPoint_00to12Z_sundy_PM10 = TotlPoint_00to12Z_sundy_file.variables['PM10'][:][:] 
TotlPoint_00to12Z_sundy_PM11 = TotlPoint_00to12Z_sundy_file.variables['PM11'][:][:] 
TotlPoint_00to12Z_sundy_PM12 = TotlPoint_00to12Z_sundy_file.variables['PM12'][:][:] 
TotlPoint_00to12Z_sundy_PM13 = TotlPoint_00to12Z_sundy_file.variables['PM13'][:][:] 
TotlPoint_00to12Z_sundy_PM14 = TotlPoint_00to12Z_sundy_file.variables['PM14'][:][:] 
TotlPoint_00to12Z_sundy_PM15 = TotlPoint_00to12Z_sundy_file.variables['PM15'][:][:] 
TotlPoint_00to12Z_sundy_PM16 = TotlPoint_00to12Z_sundy_file.variables['PM16'][:][:] 
TotlPoint_00to12Z_sundy_PM17 = TotlPoint_00to12Z_sundy_file.variables['PM17'][:][:] 
TotlPoint_00to12Z_sundy_PM18 = TotlPoint_00to12Z_sundy_file.variables['PM18'][:][:] 
TotlPoint_00to12Z_sundy_PM19 = TotlPoint_00to12Z_sundy_file.variables['PM19'][:][:] 
TotlPoint_00to12Z_sundy_Times = TotlPoint_00to12Z_sundy_file.variables['Times'][:]

###################################################################################################
#get total ROW
nROW_org, = TotlPoint_00to12Z_sundy_ITYPE.shape
nROW_extra_OG = len(LON_ng_proc)
nROW_extra = nROW_extra_OG
nROW = nROW_org + nROW_extra
print("nROW_org",nROW_org)
print("nROW_extra",nROW_extra)
print("nROW",nROW)

###################################################################################################
#Organize extra_data
extra_ITYPE_OG = np.array(ERPTYPE_ng_proc)
extra_ITYPE = extra_ITYPE_OG

extra_STKht_OG = np.array(STKHGT_ng_proc)
extra_STKht = extra_STKht_OG

extra_STKdiam_OG = np.array(STKDIAM_ng_proc)
extra_STKdiam = extra_STKdiam_OG

extra_STKtemp_OG = np.array(STKTEMP_ng_proc)
extra_STKtemp = extra_STKtemp_OG

extra_STKve_OG = np.array(STKVEL_ng_proc)
extra_STKve = extra_STKve_OG

extra_STKflw_OG = np.array(STKFLOW_ng_proc)
extra_STKflw = extra_STKflw_OG

extra_FUGht = np.empty(nROW_extra) #FUGht set as empty

extra_XLONG_OG = np.array(LON_ng_proc)
extra_XLONG = extra_XLONG_OG

extra_XLAT_OG = np.array(LAT_ng_proc)
extra_XLAT = extra_XLAT_OG

extra_STATE_OG = STATE_ng_proc

###################################################################################################
#CO2

##################################################################################
extra_CO2_FCPE_ng_proc = HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_sundy[0:12,:]
extra_CO2_OG = extra_CO2_FCPE_ng_proc

##################################################################################
extra_CO2 =extra_CO2_OG

###################################################################################################
#CH4 from OG can use GHGRP numbers

##################################################################################
extra_CH4_FCPE_ng_proc = HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_sundy[0:12,:]
extra_CH4_OG = extra_CH4_FCPE_ng_proc

###################################################################################################

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#AQ: using state-level emission ratios to CO2

###################################################################################################
#grab emission ratios from the summary ratio arrays
###################################################################################################

#based on OG process type, species, and state location
extra_X_OG = np.empty([12,nROW_extra_OG,len(species_vector)])
print("extra_X_OG.shape", extra_X_OG.shape)

for pt in range(0,nROW_extra_OG):
    proc_index = 0
    #print("proc_index",proc_index)

    state_abb_cur = extra_STATE_OG[pt]
    
    if state_abb_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_abb_cur)
        #print("state_index",state_index)
    
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * proc_spec_state_emisXdCO2_OG[proc_index,spec_index,state_index]

    else:
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * statistics.mean(proc_spec_state_emisXdCO2_OG[proc_index,spec_index,:])

extra_X_OG_dict = {}
for spec_cur in species_vector:
    #for CH4 use GHGRP numbers
    if spec_cur == 'HC01':
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_CH4_OG
    else:
        spec_index = species_vector.index(spec_cur)
        #print("spec_index",spec_index)
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_X_OG[:,:,spec_index]

###################################################################################################
#stack OG AQ species
extra_X_dict = {}
for spec_cur in species_vector:
    extra_Xi_OG = extra_X_OG_dict["extra_{0}_OG".format(spec_cur)]
    extra_X = extra_Xi_OG
    extra_X_dict["extra_{0}".format(spec_cur)] = extra_X

###################################################################################################
extra_other_spec = np.zeros([12,nROW_extra])

###################################################################################################
#append extra points to original data
TotlPoint_w_extra_00to12Z_sundy_ITYPE = np.concatenate((TotlPoint_00to12Z_sundy_ITYPE,extra_ITYPE),axis=0)
TotlPoint_w_extra_00to12Z_sundy_STKht = np.concatenate((TotlPoint_00to12Z_sundy_STKht,extra_STKht),axis=0)
TotlPoint_w_extra_00to12Z_sundy_STKdiam = np.concatenate((TotlPoint_00to12Z_sundy_STKdiam,extra_STKdiam),axis=0)
TotlPoint_w_extra_00to12Z_sundy_STKtemp = np.concatenate((TotlPoint_00to12Z_sundy_STKtemp,extra_STKtemp),axis=0)
TotlPoint_w_extra_00to12Z_sundy_STKve = np.concatenate((TotlPoint_00to12Z_sundy_STKve,extra_STKve),axis=0)
TotlPoint_w_extra_00to12Z_sundy_STKflw = np.concatenate((TotlPoint_00to12Z_sundy_STKflw,extra_STKflw),axis=0)
TotlPoint_w_extra_00to12Z_sundy_FUGht = np.concatenate((TotlPoint_00to12Z_sundy_FUGht,extra_FUGht),axis=0)
TotlPoint_w_extra_00to12Z_sundy_XLONG = np.concatenate((TotlPoint_00to12Z_sundy_XLONG,extra_XLONG),axis=0)
TotlPoint_w_extra_00to12Z_sundy_XLAT = np.concatenate((TotlPoint_00to12Z_sundy_XLAT,extra_XLAT),axis=0)
TotlPoint_w_extra_00to12Z_sundy_CO2 = np.concatenate((TotlPoint_00to12Z_sundy_CO2,extra_CO2),axis=1)
TotlPoint_w_extra_00to12Z_sundy_CO = np.concatenate((TotlPoint_00to12Z_sundy_CO,extra_X_dict["extra_CO"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_NH3 = np.concatenate((TotlPoint_00to12Z_sundy_NH3,extra_X_dict["extra_NH3"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_NOX = np.concatenate((TotlPoint_00to12Z_sundy_NOX,extra_X_dict["extra_NOX"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM10_PRI = np.concatenate((TotlPoint_00to12Z_sundy_PM10_PRI,extra_X_dict["extra_PM10-PRI"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM25_PRI = np.concatenate((TotlPoint_00to12Z_sundy_PM25_PRI,extra_X_dict["extra_PM25-PRI"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_SO2 = np.concatenate((TotlPoint_00to12Z_sundy_SO2,extra_X_dict["extra_SO2"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_VOC = np.concatenate((TotlPoint_00to12Z_sundy_VOC,extra_X_dict["extra_VOC"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC01 = np.concatenate((TotlPoint_00to12Z_sundy_HC01,extra_X_dict["extra_HC01"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC02 = np.concatenate((TotlPoint_00to12Z_sundy_HC02,extra_X_dict["extra_HC02"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC03 = np.concatenate((TotlPoint_00to12Z_sundy_HC03,extra_X_dict["extra_HC03"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC04 = np.concatenate((TotlPoint_00to12Z_sundy_HC04,extra_X_dict["extra_HC04"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC05 = np.concatenate((TotlPoint_00to12Z_sundy_HC05,extra_X_dict["extra_HC05"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC06 = np.concatenate((TotlPoint_00to12Z_sundy_HC06,extra_X_dict["extra_HC06"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC07 = np.concatenate((TotlPoint_00to12Z_sundy_HC07,extra_X_dict["extra_HC07"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC08 = np.concatenate((TotlPoint_00to12Z_sundy_HC08,extra_X_dict["extra_HC08"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC09 = np.concatenate((TotlPoint_00to12Z_sundy_HC09,extra_X_dict["extra_HC09"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC10 = np.concatenate((TotlPoint_00to12Z_sundy_HC10,extra_X_dict["extra_HC10"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC11 = np.concatenate((TotlPoint_00to12Z_sundy_HC11,extra_X_dict["extra_HC11"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC12 = np.concatenate((TotlPoint_00to12Z_sundy_HC12,extra_X_dict["extra_HC12"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC13 = np.concatenate((TotlPoint_00to12Z_sundy_HC13,extra_X_dict["extra_HC13"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC14 = np.concatenate((TotlPoint_00to12Z_sundy_HC14,extra_X_dict["extra_HC14"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC15 = np.concatenate((TotlPoint_00to12Z_sundy_HC15,extra_X_dict["extra_HC15"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC16 = np.concatenate((TotlPoint_00to12Z_sundy_HC16,extra_X_dict["extra_HC16"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC17 = np.concatenate((TotlPoint_00to12Z_sundy_HC17,extra_X_dict["extra_HC17"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC18 = np.concatenate((TotlPoint_00to12Z_sundy_HC18,extra_X_dict["extra_HC18"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC19 = np.concatenate((TotlPoint_00to12Z_sundy_HC19,extra_X_dict["extra_HC19"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC20 = np.concatenate((TotlPoint_00to12Z_sundy_HC20,extra_X_dict["extra_HC20"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC21 = np.concatenate((TotlPoint_00to12Z_sundy_HC21,extra_X_dict["extra_HC21"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC22 = np.concatenate((TotlPoint_00to12Z_sundy_HC22,extra_X_dict["extra_HC22"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC23 = np.concatenate((TotlPoint_00to12Z_sundy_HC23,extra_X_dict["extra_HC23"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC24 = np.concatenate((TotlPoint_00to12Z_sundy_HC24,extra_X_dict["extra_HC24"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC25 = np.concatenate((TotlPoint_00to12Z_sundy_HC25,extra_X_dict["extra_HC25"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC26 = np.concatenate((TotlPoint_00to12Z_sundy_HC26,extra_X_dict["extra_HC26"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC27 = np.concatenate((TotlPoint_00to12Z_sundy_HC27,extra_X_dict["extra_HC27"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC28 = np.concatenate((TotlPoint_00to12Z_sundy_HC28,extra_X_dict["extra_HC28"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC29 = np.concatenate((TotlPoint_00to12Z_sundy_HC29,extra_X_dict["extra_HC29"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC30 = np.concatenate((TotlPoint_00to12Z_sundy_HC30,extra_X_dict["extra_HC30"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC31 = np.concatenate((TotlPoint_00to12Z_sundy_HC31,extra_X_dict["extra_HC31"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC32 = np.concatenate((TotlPoint_00to12Z_sundy_HC32,extra_X_dict["extra_HC32"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC33 = np.concatenate((TotlPoint_00to12Z_sundy_HC33,extra_X_dict["extra_HC33"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC34 = np.concatenate((TotlPoint_00to12Z_sundy_HC34,extra_X_dict["extra_HC34"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC35 = np.concatenate((TotlPoint_00to12Z_sundy_HC35,extra_X_dict["extra_HC35"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC36 = np.concatenate((TotlPoint_00to12Z_sundy_HC36,extra_X_dict["extra_HC36"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC37 = np.concatenate((TotlPoint_00to12Z_sundy_HC37,extra_X_dict["extra_HC37"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC38 = np.concatenate((TotlPoint_00to12Z_sundy_HC38,extra_X_dict["extra_HC38"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC39 = np.concatenate((TotlPoint_00to12Z_sundy_HC39,extra_X_dict["extra_HC39"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC40 = np.concatenate((TotlPoint_00to12Z_sundy_HC40,extra_X_dict["extra_HC40"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC41 = np.concatenate((TotlPoint_00to12Z_sundy_HC41,extra_X_dict["extra_HC41"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC42 = np.concatenate((TotlPoint_00to12Z_sundy_HC42,extra_X_dict["extra_HC42"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC43 = np.concatenate((TotlPoint_00to12Z_sundy_HC43,extra_X_dict["extra_HC43"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC44 = np.concatenate((TotlPoint_00to12Z_sundy_HC44,extra_X_dict["extra_HC44"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC45 = np.concatenate((TotlPoint_00to12Z_sundy_HC45,extra_X_dict["extra_HC45"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC46 = np.concatenate((TotlPoint_00to12Z_sundy_HC46,extra_X_dict["extra_HC46"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC47 = np.concatenate((TotlPoint_00to12Z_sundy_HC47,extra_X_dict["extra_HC47"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC48 = np.concatenate((TotlPoint_00to12Z_sundy_HC48,extra_X_dict["extra_HC48"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC49 = np.concatenate((TotlPoint_00to12Z_sundy_HC49,extra_X_dict["extra_HC49"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC50 = np.concatenate((TotlPoint_00to12Z_sundy_HC50,extra_X_dict["extra_HC50"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC51 = np.concatenate((TotlPoint_00to12Z_sundy_HC51,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC52 = np.concatenate((TotlPoint_00to12Z_sundy_HC52,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC53 = np.concatenate((TotlPoint_00to12Z_sundy_HC53,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC54 = np.concatenate((TotlPoint_00to12Z_sundy_HC54,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC55 = np.concatenate((TotlPoint_00to12Z_sundy_HC55,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC56 = np.concatenate((TotlPoint_00to12Z_sundy_HC56,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC57 = np.concatenate((TotlPoint_00to12Z_sundy_HC57,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC58 = np.concatenate((TotlPoint_00to12Z_sundy_HC58,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC59 = np.concatenate((TotlPoint_00to12Z_sundy_HC59,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC60 = np.concatenate((TotlPoint_00to12Z_sundy_HC60,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC61 = np.concatenate((TotlPoint_00to12Z_sundy_HC61,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC62 = np.concatenate((TotlPoint_00to12Z_sundy_HC62,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC63 = np.concatenate((TotlPoint_00to12Z_sundy_HC63,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC64 = np.concatenate((TotlPoint_00to12Z_sundy_HC64,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC65 = np.concatenate((TotlPoint_00to12Z_sundy_HC65,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC66 = np.concatenate((TotlPoint_00to12Z_sundy_HC66,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC67 = np.concatenate((TotlPoint_00to12Z_sundy_HC67,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_HC68 = np.concatenate((TotlPoint_00to12Z_sundy_HC68,extra_other_spec),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM01 = np.concatenate((TotlPoint_00to12Z_sundy_PM01,extra_X_dict["extra_PM01"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM02 = np.concatenate((TotlPoint_00to12Z_sundy_PM02,extra_X_dict["extra_PM02"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM03 = np.concatenate((TotlPoint_00to12Z_sundy_PM03,extra_X_dict["extra_PM03"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM04 = np.concatenate((TotlPoint_00to12Z_sundy_PM04,extra_X_dict["extra_PM04"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM05 = np.concatenate((TotlPoint_00to12Z_sundy_PM05,extra_X_dict["extra_PM05"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM06 = np.concatenate((TotlPoint_00to12Z_sundy_PM06,extra_X_dict["extra_PM06"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM07 = np.concatenate((TotlPoint_00to12Z_sundy_PM07,extra_X_dict["extra_PM07"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM08 = np.concatenate((TotlPoint_00to12Z_sundy_PM08,extra_X_dict["extra_PM08"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM09 = np.concatenate((TotlPoint_00to12Z_sundy_PM09,extra_X_dict["extra_PM09"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM10 = np.concatenate((TotlPoint_00to12Z_sundy_PM10,extra_X_dict["extra_PM10"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM11 = np.concatenate((TotlPoint_00to12Z_sundy_PM11,extra_X_dict["extra_PM11"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM12 = np.concatenate((TotlPoint_00to12Z_sundy_PM12,extra_X_dict["extra_PM12"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM13 = np.concatenate((TotlPoint_00to12Z_sundy_PM13,extra_X_dict["extra_PM13"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM14 = np.concatenate((TotlPoint_00to12Z_sundy_PM14,extra_X_dict["extra_PM14"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM15 = np.concatenate((TotlPoint_00to12Z_sundy_PM15,extra_X_dict["extra_PM15"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM16 = np.concatenate((TotlPoint_00to12Z_sundy_PM16,extra_X_dict["extra_PM16"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM17 = np.concatenate((TotlPoint_00to12Z_sundy_PM17,extra_X_dict["extra_PM17"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM18 = np.concatenate((TotlPoint_00to12Z_sundy_PM18,extra_X_dict["extra_PM18"]),axis=1)
TotlPoint_w_extra_00to12Z_sundy_PM19 = np.concatenate((TotlPoint_00to12Z_sundy_PM19,extra_X_dict["extra_PM19"]),axis=1)

###################################################################################################
#write total points with extra points appended
TotlPoint_w_extra_00to12Z_sundy_fn = append_dir+'/sundy/PtOG_00to12Z.nc'
TotlPoint_w_extra_00to12Z_sundy_file = Dataset(TotlPoint_w_extra_00to12Z_sundy_fn,mode='w',format='NETCDF3_64BIT')

#Creat dimensions
TotlPoint_w_extra_00to12Z_sundy_file.createDimension("ROW", nROW)
TotlPoint_w_extra_00to12Z_sundy_file.createDimension("Time", 12)
TotlPoint_w_extra_00to12Z_sundy_file.sync()

#Create variables
#float ITYPE(ROW) ;
ITYPE = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('ITYPE','f4',('ROW'),fill_value = float(0))
ITYPE[:] = TotlPoint_w_extra_00to12Z_sundy_ITYPE
varattrs=["FieldType","MemoryOrder","description","units","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['ITYPE'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['ITYPE'], varattr);
        setattr(ITYPE, varattr, varattrVal)

#float STKht(ROW) ;
STKht = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('STKht','f4',('ROW'),fill_value = float(0))
STKht[:] = TotlPoint_w_extra_00to12Z_sundy_STKht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['STKht'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['STKht'], varattr);
        setattr(STKht, varattr, varattrVal)
        
#float STKdiam(ROW) ;
STKdiam = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('STKdiam','f4',('ROW'),fill_value = float(0))
STKdiam[:] = TotlPoint_w_extra_00to12Z_sundy_STKdiam
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['STKdiam'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['STKdiam'], varattr);
        setattr(STKdiam, varattr, varattrVal)

#float STKtemp(ROW) ;
STKtemp = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('STKtemp','f4',('ROW'),fill_value = float(0))
STKtemp[:] = TotlPoint_w_extra_00to12Z_sundy_STKtemp
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['STKtemp'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['STKtemp'], varattr);
        setattr(STKtemp, varattr, varattrVal)
        
#float STKve(ROW) ;
STKve = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('STKve','f4',('ROW'),fill_value = float(0))
STKve[:] = TotlPoint_w_extra_00to12Z_sundy_STKve
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['STKve'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['STKve'], varattr);
        setattr(STKve, varattr, varattrVal)
        
#float STKflw(ROW) ;
STKflw = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('STKflw','f4',('ROW'),fill_value = float(0))
STKflw[:] = TotlPoint_w_extra_00to12Z_sundy_STKflw
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['STKflw'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['STKflw'], varattr);
        setattr(STKflw, varattr, varattrVal)
        
#float FUGht(ROW) ;
FUGht = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('FUGht','f4',('ROW'),fill_value = float(0))
FUGht[:] = TotlPoint_w_extra_00to12Z_sundy_FUGht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['FUGht'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['FUGht'], varattr);
        setattr(FUGht, varattr, varattrVal)
        
#float XLONG(ROW) ;
XLONG = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('XLONG','f4',('ROW'),fill_value = float(0))
XLONG[:] = TotlPoint_w_extra_00to12Z_sundy_XLONG
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['XLONG'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['XLONG'], varattr);
        setattr(XLONG, varattr, varattrVal)
        
#float XLAT(ROW) ;
XLAT = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('XLAT','f4',('ROW'),fill_value = float(0))
XLAT[:] = TotlPoint_w_extra_00to12Z_sundy_XLAT
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['XLAT'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['XLAT'], varattr);
        setattr(XLAT, varattr, varattrVal)

#float CO2(Time, ROW) ;
CO2 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('CO2','f4',('Time','ROW'))
CO2[:] = TotlPoint_w_extra_00to12Z_sundy_CO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['CO2'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['CO2'], varattr);
        setattr(CO2, varattr, varattrVal)
        
#float CO(Time, ROW) ;
CO = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('CO','f4',('Time','ROW'))
CO[:] = TotlPoint_w_extra_00to12Z_sundy_CO
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['CO'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['CO'], varattr);
        setattr(CO, varattr, varattrVal)
        
#float NH3(Time, ROW) ;
NH3 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('NH3','f4',('Time','ROW'))
NH3[:] = TotlPoint_w_extra_00to12Z_sundy_NH3
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['NH3'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['NH3'], varattr);
        setattr(NH3, varattr, varattrVal)
        
#float NOX(Time, ROW) ;
NOX = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('NOX','f4',('Time','ROW'))
NOX[:] = TotlPoint_w_extra_00to12Z_sundy_NOX
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['NOX'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['NOX'], varattr);
        setattr(NOX, varattr, varattrVal)
        
#float PM10-PRI(Time, ROW) ;
PM10_PRI = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM10-PRI','f4',('Time','ROW'))
PM10_PRI[:] = TotlPoint_w_extra_00to12Z_sundy_PM10_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM10-PRI'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM10-PRI'], varattr);
        setattr(PM10_PRI, varattr, varattrVal)
        
#float PM25-PRI(Time, ROW) ;
PM25_PRI = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM25-PRI','f4',('Time','ROW'))
PM25_PRI[:] = TotlPoint_w_extra_00to12Z_sundy_PM25_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM25-PRI'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM25-PRI'], varattr);
        setattr(PM25_PRI, varattr, varattrVal)
        
#float SO2(Time, ROW) ;
SO2 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('SO2','f4',('Time','ROW'))
SO2[:] = TotlPoint_w_extra_00to12Z_sundy_SO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['SO2'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['SO2'], varattr);
        setattr(SO2, varattr, varattrVal)
        
#float VOC(Time, ROW) ;
VOC = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('VOC','f4',('Time','ROW'))
VOC[:] = TotlPoint_w_extra_00to12Z_sundy_VOC
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['VOC'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['VOC'], varattr);
        setattr(VOC, varattr, varattrVal)
        
#float HC01(Time, ROW) ;
HC01 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC01','f4',('Time','ROW'))
HC01[:] = TotlPoint_w_extra_00to12Z_sundy_HC01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC01'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC01'], varattr);
        setattr(HC01, varattr, varattrVal)
        
#float HC02(Time, ROW) ;
HC02 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC02','f4',('Time','ROW'))
HC02[:] = TotlPoint_w_extra_00to12Z_sundy_HC02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC02'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC02'], varattr);
        setattr(HC02, varattr, varattrVal)
        
#float HC03(Time, ROW) ;
HC03 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC03','f4',('Time','ROW'))
HC03[:] = TotlPoint_w_extra_00to12Z_sundy_HC03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC03'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC03'], varattr);
        setattr(HC03, varattr, varattrVal)
        
#float HC04(Time, ROW) ;
HC04 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC04','f4',('Time','ROW'))
HC04[:] = TotlPoint_w_extra_00to12Z_sundy_HC04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC04'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC04'], varattr);
        setattr(HC04, varattr, varattrVal)
        
#float HC05(Time, ROW) ;
HC05 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC05','f4',('Time','ROW'))
HC05[:] = TotlPoint_w_extra_00to12Z_sundy_HC05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC05'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC05'], varattr);
        setattr(HC05, varattr, varattrVal)
        
#float HC06(Time, ROW) ;
HC06 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC06','f4',('Time','ROW'))
HC06[:] = TotlPoint_w_extra_00to12Z_sundy_HC06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC06'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC06'], varattr);
        setattr(HC06, varattr, varattrVal)
        
#float HC07(Time, ROW) ;
HC07 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC07','f4',('Time','ROW'))
HC07[:] = TotlPoint_w_extra_00to12Z_sundy_HC07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC07'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC07'], varattr);
        setattr(HC07, varattr, varattrVal)
        
#float HC08(Time, ROW) ;
HC08 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC08','f4',('Time','ROW'))
HC08[:] = TotlPoint_w_extra_00to12Z_sundy_HC08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC08'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC08'], varattr);
        setattr(HC08, varattr, varattrVal)
        
#float HC09(Time, ROW) ;
HC09 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC09','f4',('Time','ROW'))
HC09[:] = TotlPoint_w_extra_00to12Z_sundy_HC09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC09'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC09'], varattr);
        setattr(HC09, varattr, varattrVal)
        
#float HC10(Time, ROW) ;
HC10 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC10','f4',('Time','ROW'))
HC10[:] = TotlPoint_w_extra_00to12Z_sundy_HC10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC10'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC10'], varattr);
        setattr(HC10, varattr, varattrVal)

#float HC11(Time, ROW) ;
HC11 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC11','f4',('Time','ROW'))
HC11[:] = TotlPoint_w_extra_00to12Z_sundy_HC11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC11'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC11'], varattr);
        setattr(HC11, varattr, varattrVal)
        
#float HC12(Time, ROW) ;
HC12 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC12','f4',('Time','ROW'))
HC12[:] = TotlPoint_w_extra_00to12Z_sundy_HC12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC12'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC12'], varattr);
        setattr(HC12, varattr, varattrVal)
        
#float HC13(Time, ROW) ;
HC13 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC13','f4',('Time','ROW'))
HC13[:] = TotlPoint_w_extra_00to12Z_sundy_HC13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC13'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC13'], varattr);
        setattr(HC13, varattr, varattrVal)
        
#float HC14(Time, ROW) ;
HC14 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC14','f4',('Time','ROW'))
HC14[:] = TotlPoint_w_extra_00to12Z_sundy_HC14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC14'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC14'], varattr);
        setattr(HC14, varattr, varattrVal)
        
#float HC15(Time, ROW) ;
HC15 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC15','f4',('Time','ROW'))
HC15[:] = TotlPoint_w_extra_00to12Z_sundy_HC15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC15'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC15'], varattr);
        setattr(HC15, varattr, varattrVal)
        
#float HC16(Time, ROW) ;
HC16 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC16','f4',('Time','ROW'))
HC16[:] = TotlPoint_w_extra_00to12Z_sundy_HC16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC16'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC16'], varattr);
        setattr(HC16, varattr, varattrVal)
        
#float HC17(Time, ROW) ;
HC17 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC17','f4',('Time','ROW'))
HC17[:] = TotlPoint_w_extra_00to12Z_sundy_HC17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC17'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC17'], varattr);
        setattr(HC17, varattr, varattrVal)
        
#float HC18(Time, ROW) ;
HC18 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC18','f4',('Time','ROW'))
HC18[:] = TotlPoint_w_extra_00to12Z_sundy_HC18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC18'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC18'], varattr);
        setattr(HC18, varattr, varattrVal)
        
#float HC19(Time, ROW) ;
HC19 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC19','f4',('Time','ROW'))
HC19[:] = TotlPoint_w_extra_00to12Z_sundy_HC19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC19'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC19'], varattr);
        setattr(HC19, varattr, varattrVal)

#float HC20(Time, ROW) ;
HC20 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC20','f4',('Time','ROW'))
HC20[:] = TotlPoint_w_extra_00to12Z_sundy_HC20
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC20'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC20'], varattr);
        setattr(HC20, varattr, varattrVal)

#float HC21(Time, ROW) ;
HC21 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC21','f4',('Time','ROW'))
HC21[:] = TotlPoint_w_extra_00to12Z_sundy_HC21
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC21'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC21'], varattr);
        setattr(HC21, varattr, varattrVal)
        
#float HC22(Time, ROW) ;
HC22 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC22','f4',('Time','ROW'))
HC22[:] = TotlPoint_w_extra_00to12Z_sundy_HC22
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC22'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC22'], varattr);
        setattr(HC22, varattr, varattrVal)
        
#float HC23(Time, ROW) ;
HC23 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC23','f4',('Time','ROW'))
HC23[:] = TotlPoint_w_extra_00to12Z_sundy_HC23
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC23'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC23'], varattr);
        setattr(HC23, varattr, varattrVal)
        
#float HC24(Time, ROW) ;
HC24 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC24','f4',('Time','ROW'))
HC24[:] = TotlPoint_w_extra_00to12Z_sundy_HC24
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC24'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC24'], varattr);
        setattr(HC24, varattr, varattrVal)
        
#float HC25(Time, ROW) ;
HC25 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC25','f4',('Time','ROW'))
HC25[:] = TotlPoint_w_extra_00to12Z_sundy_HC25
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC25'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC25'], varattr);
        setattr(HC25, varattr, varattrVal)
        
#float HC26(Time, ROW) ;
HC26 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC26','f4',('Time','ROW'))
HC26[:] = TotlPoint_w_extra_00to12Z_sundy_HC26
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC26'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC26'], varattr);
        setattr(HC26, varattr, varattrVal)
        
#float HC27(Time, ROW) ;
HC27 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC27','f4',('Time','ROW'))
HC27[:] = TotlPoint_w_extra_00to12Z_sundy_HC27
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC27'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC27'], varattr);
        setattr(HC27, varattr, varattrVal)
        
#float HC28(Time, ROW) ;
HC28 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC28','f4',('Time','ROW'))
HC28[:] = TotlPoint_w_extra_00to12Z_sundy_HC28
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC28'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC28'], varattr);
        setattr(HC28, varattr, varattrVal)
        
#float HC29(Time, ROW) ;
HC29 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC29','f4',('Time','ROW'))
HC29[:] = TotlPoint_w_extra_00to12Z_sundy_HC29
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC29'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC29'], varattr);
        setattr(HC29, varattr, varattrVal)

#float HC30(Time, ROW) ;
HC30 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC30','f4',('Time','ROW'))
HC30[:] = TotlPoint_w_extra_00to12Z_sundy_HC30
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC30'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC30'], varattr);
        setattr(HC30, varattr, varattrVal)

#float HC31(Time, ROW) ;
HC31 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC31','f4',('Time','ROW'))
HC31[:] = TotlPoint_w_extra_00to12Z_sundy_HC31
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC31'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC31'], varattr);
        setattr(HC31, varattr, varattrVal)
        
#float HC32(Time, ROW) ;
HC32 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC32','f4',('Time','ROW'))
HC32[:] = TotlPoint_w_extra_00to12Z_sundy_HC32
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC32'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC32'], varattr);
        setattr(HC32, varattr, varattrVal)
        
#float HC33(Time, ROW) ;
HC33 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC33','f4',('Time','ROW'))
HC33[:] = TotlPoint_w_extra_00to12Z_sundy_HC33
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC33'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC33'], varattr);
        setattr(HC33, varattr, varattrVal)
        
#float HC34(Time, ROW) ;
HC34 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC34','f4',('Time','ROW'))
HC34[:] = TotlPoint_w_extra_00to12Z_sundy_HC34
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC34'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC34'], varattr);
        setattr(HC34, varattr, varattrVal)
        
#float HC35(Time, ROW) ;
HC35 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC35','f4',('Time','ROW'))
HC35[:] = TotlPoint_w_extra_00to12Z_sundy_HC35
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC35'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC35'], varattr);
        setattr(HC35, varattr, varattrVal)
        
#float HC36(Time, ROW) ;
HC36 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC36','f4',('Time','ROW'))
HC36[:] = TotlPoint_w_extra_00to12Z_sundy_HC36
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC36'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC36'], varattr);
        setattr(HC36, varattr, varattrVal)
        
#float HC37(Time, ROW) ;
HC37 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC37','f4',('Time','ROW'))
HC37[:] = TotlPoint_w_extra_00to12Z_sundy_HC37
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC37'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC37'], varattr);
        setattr(HC37, varattr, varattrVal)
        
#float HC38(Time, ROW) ;
HC38 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC38','f4',('Time','ROW'))
HC38[:] = TotlPoint_w_extra_00to12Z_sundy_HC38
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC38'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC38'], varattr);
        setattr(HC38, varattr, varattrVal)
        
#float HC39(Time, ROW) ;
HC39 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC39','f4',('Time','ROW'))
HC39[:] = TotlPoint_w_extra_00to12Z_sundy_HC39
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC39'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC39'], varattr);
        setattr(HC39, varattr, varattrVal)
        
#float HC40(Time, ROW) ;
HC40 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC40','f4',('Time','ROW'))
HC40[:] = TotlPoint_w_extra_00to12Z_sundy_HC40
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC40'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC40'], varattr);
        setattr(HC40, varattr, varattrVal)

#float HC41(Time, ROW) ;
HC41 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC41','f4',('Time','ROW'))
HC41[:] = TotlPoint_w_extra_00to12Z_sundy_HC41
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC41'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC41'], varattr);
        setattr(HC41, varattr, varattrVal)
        
#float HC42(Time, ROW) ;
HC42 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC42','f4',('Time','ROW'))
HC42[:] = TotlPoint_w_extra_00to12Z_sundy_HC42
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC42'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC42'], varattr);
        setattr(HC42, varattr, varattrVal)
        
#float HC43(Time, ROW) ;
HC43 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC43','f4',('Time','ROW'))
HC43[:] = TotlPoint_w_extra_00to12Z_sundy_HC43
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC43'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC43'], varattr);
        setattr(HC43, varattr, varattrVal)
        
#float HC44(Time, ROW) ;
HC44 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC44','f4',('Time','ROW'))
HC44[:] = TotlPoint_w_extra_00to12Z_sundy_HC44
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC44'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC44'], varattr);
        setattr(HC44, varattr, varattrVal)
        
#float HC45(Time, ROW) ;
HC45 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC45','f4',('Time','ROW'))
HC45[:] = TotlPoint_w_extra_00to12Z_sundy_HC45
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC45'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC45'], varattr);
        setattr(HC45, varattr, varattrVal)
        
#float HC46(Time, ROW) ;
HC46 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC46','f4',('Time','ROW'))
HC46[:] = TotlPoint_w_extra_00to12Z_sundy_HC46
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC46'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC46'], varattr);
        setattr(HC46, varattr, varattrVal)
        
#float HC47(Time, ROW) ;
HC47 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC47','f4',('Time','ROW'))
HC47[:] = TotlPoint_w_extra_00to12Z_sundy_HC47
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC47'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC47'], varattr);
        setattr(HC47, varattr, varattrVal)
        
#float HC48(Time, ROW) ;
HC48 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC48','f4',('Time','ROW'))
HC48[:] = TotlPoint_w_extra_00to12Z_sundy_HC48
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC48'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC48'], varattr);
        setattr(HC48, varattr, varattrVal)
        
#float HC49(Time, ROW) ;
HC49 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC49','f4',('Time','ROW'))
HC49[:] = TotlPoint_w_extra_00to12Z_sundy_HC49
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC49'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC49'], varattr);
        setattr(HC49, varattr, varattrVal)
        
#float HC50(Time, ROW) ;
HC50 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC50','f4',('Time','ROW'))
HC50[:] = TotlPoint_w_extra_00to12Z_sundy_HC50
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC50'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC50'], varattr);
        setattr(HC50, varattr, varattrVal)

#float HC51(Time, ROW) ;
HC51 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC51','f4',('Time','ROW'))
HC51[:] = TotlPoint_w_extra_00to12Z_sundy_HC51
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC51'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC51'], varattr);
        setattr(HC51, varattr, varattrVal)
        
#float HC52(Time, ROW) ;
HC52 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC52','f4',('Time','ROW'))
HC52[:] = TotlPoint_w_extra_00to12Z_sundy_HC52
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC52'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC52'], varattr);
        setattr(HC52, varattr, varattrVal)
        
#float HC53(Time, ROW) ;
HC53 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC53','f4',('Time','ROW'))
HC53[:] = TotlPoint_w_extra_00to12Z_sundy_HC53
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC53'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC53'], varattr);
        setattr(HC53, varattr, varattrVal)
        
#float HC54(Time, ROW) ;
HC54 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC54','f4',('Time','ROW'))
HC54[:] = TotlPoint_w_extra_00to12Z_sundy_HC54
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC54'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC54'], varattr);
        setattr(HC54, varattr, varattrVal)
        
#float HC55(Time, ROW) ;
HC55 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC55','f4',('Time','ROW'))
HC55[:] = TotlPoint_w_extra_00to12Z_sundy_HC55
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC55'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC55'], varattr);
        setattr(HC55, varattr, varattrVal)
        
#float HC56(Time, ROW) ;
HC56 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC56','f4',('Time','ROW'))
HC56[:] = TotlPoint_w_extra_00to12Z_sundy_HC56
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC56'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC56'], varattr);
        setattr(HC56, varattr, varattrVal)
        
#float HC57(Time, ROW) ;
HC57 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC57','f4',('Time','ROW'))
HC57[:] = TotlPoint_w_extra_00to12Z_sundy_HC57
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC57'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC57'], varattr);
        setattr(HC57, varattr, varattrVal)
        
#float HC58(Time, ROW) ;
HC58 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC58','f4',('Time','ROW'))
HC58[:] = TotlPoint_w_extra_00to12Z_sundy_HC58
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC58'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC58'], varattr);
        setattr(HC58, varattr, varattrVal)
        
#float HC59(Time, ROW) ;
HC59 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC59','f4',('Time','ROW'))
HC59[:] = TotlPoint_w_extra_00to12Z_sundy_HC59
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC59'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC59'], varattr);
        setattr(HC59, varattr, varattrVal)

#float HC60(Time, ROW) ;
HC60 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC60','f4',('Time','ROW'))
HC60[:] = TotlPoint_w_extra_00to12Z_sundy_HC60
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC60'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC60'], varattr);
        setattr(HC60, varattr, varattrVal)

#float HC61(Time, ROW) ;
HC61 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC61','f4',('Time','ROW'))
HC61[:] = TotlPoint_w_extra_00to12Z_sundy_HC61
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC61'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC61'], varattr);
        setattr(HC61, varattr, varattrVal)
        
#float HC62(Time, ROW) ;
HC62 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC62','f4',('Time','ROW'))
HC62[:] = TotlPoint_w_extra_00to12Z_sundy_HC62
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC62'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC62'], varattr);
        setattr(HC62, varattr, varattrVal)
        
#float HC63(Time, ROW) ;
HC63 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC63','f4',('Time','ROW'))
HC63[:] = TotlPoint_w_extra_00to12Z_sundy_HC63
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC63'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC63'], varattr);
        setattr(HC63, varattr, varattrVal)
        
#float HC64(Time, ROW) ;
HC64 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC64','f4',('Time','ROW'))
HC64[:] = TotlPoint_w_extra_00to12Z_sundy_HC64
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC64'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC64'], varattr);
        setattr(HC64, varattr, varattrVal)
        
#float HC65(Time, ROW) ;
HC65 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC65','f4',('Time','ROW'))
HC65[:] = TotlPoint_w_extra_00to12Z_sundy_HC65
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC65'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC65'], varattr);
        setattr(HC65, varattr, varattrVal)
        
#float HC66(Time, ROW) ;
HC66 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC66','f4',('Time','ROW'))
HC66[:] = TotlPoint_w_extra_00to12Z_sundy_HC66
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC66'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC66'], varattr);
        setattr(HC66, varattr, varattrVal)
        
#float HC67(Time, ROW) ;
HC67 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC67','f4',('Time','ROW'))
HC67[:] = TotlPoint_w_extra_00to12Z_sundy_HC67
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC67'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC67'], varattr);
        setattr(HC67, varattr, varattrVal)
        
#float HC68(Time, ROW) ;
HC68 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('HC68','f4',('Time','ROW'))
HC68[:] = TotlPoint_w_extra_00to12Z_sundy_HC68
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['HC68'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['HC68'], varattr);
        setattr(HC68, varattr, varattrVal)

#float PM01(Time, ROW) ;
PM01 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM01','f4',('Time','ROW'))
PM01[:] = TotlPoint_w_extra_00to12Z_sundy_PM01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM01'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM01'], varattr);
        setattr(PM01, varattr, varattrVal)

#float PM02(Time, ROW) ;
PM02 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM02','f4',('Time','ROW'))
PM02[:] = TotlPoint_w_extra_00to12Z_sundy_PM02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM02'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM02'], varattr);
        setattr(PM02, varattr, varattrVal)
        
#float PM03(Time, ROW) ;
PM03 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM03','f4',('Time','ROW'))
PM03[:] = TotlPoint_w_extra_00to12Z_sundy_PM03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM03'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM03'], varattr);
        setattr(PM03, varattr, varattrVal)

#float PM04(Time, ROW) ;
PM04 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM04','f4',('Time','ROW'))
PM04[:] = TotlPoint_w_extra_00to12Z_sundy_PM04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM04'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM04'], varattr);
        setattr(PM04, varattr, varattrVal)
        
#float PM05(Time, ROW) ;
PM05 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM05','f4',('Time','ROW'))
PM05[:] = TotlPoint_w_extra_00to12Z_sundy_PM05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM05'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM05'], varattr);
        setattr(PM05, varattr, varattrVal)

#float PM06(Time, ROW) ;
PM06 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM06','f4',('Time','ROW'))
PM06[:] = TotlPoint_w_extra_00to12Z_sundy_PM06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM06'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM06'], varattr);
        setattr(PM06, varattr, varattrVal)
        
#float PM07(Time, ROW) ;
PM07 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM07','f4',('Time','ROW'))
PM07[:] = TotlPoint_w_extra_00to12Z_sundy_PM07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM07'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM07'], varattr);
        setattr(PM07, varattr, varattrVal)
        
#float PM08(Time, ROW) ;
PM08 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM08','f4',('Time','ROW'))
PM08[:] = TotlPoint_w_extra_00to12Z_sundy_PM08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM08'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM08'], varattr);
        setattr(PM08, varattr, varattrVal)
        
#float PM09(Time, ROW) ;
PM09 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM09','f4',('Time','ROW'))
PM09[:] = TotlPoint_w_extra_00to12Z_sundy_PM09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM09'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM09'], varattr);
        setattr(PM09, varattr, varattrVal)
        
#float PM10(Time, ROW) ;
PM10 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM10','f4',('Time','ROW'))
PM10[:] = TotlPoint_w_extra_00to12Z_sundy_PM10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM10'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM10'], varattr);
        setattr(PM10, varattr, varattrVal)
        
#float PM11(Time, ROW) ;
PM11 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM11','f4',('Time','ROW'))
PM11[:] = TotlPoint_w_extra_00to12Z_sundy_PM11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM11'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM11'], varattr);
        setattr(PM11, varattr, varattrVal)

#float PM12(Time, ROW) ;
PM12 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM12','f4',('Time','ROW'))
PM12[:] = TotlPoint_w_extra_00to12Z_sundy_PM12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM12'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM12'], varattr);
        setattr(PM12, varattr, varattrVal)
        
#float PM13(Time, ROW) ;
PM13 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM13','f4',('Time','ROW'))
PM13[:] = TotlPoint_w_extra_00to12Z_sundy_PM13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM13'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM13'], varattr);
        setattr(PM13, varattr, varattrVal)

#float PM14(Time, ROW) ;
PM14 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM14','f4',('Time','ROW'))
PM14[:] = TotlPoint_w_extra_00to12Z_sundy_PM14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM14'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM14'], varattr);
        setattr(PM14, varattr, varattrVal)
        
#float PM15(Time, ROW) ;
PM15 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM15','f4',('Time','ROW'))
PM15[:] = TotlPoint_w_extra_00to12Z_sundy_PM15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM15'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM15'], varattr);
        setattr(PM15, varattr, varattrVal)

#float PM16(Time, ROW) ;
PM16 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM16','f4',('Time','ROW'))
PM16[:] = TotlPoint_w_extra_00to12Z_sundy_PM16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM16'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM16'], varattr);
        setattr(PM16, varattr, varattrVal)
        
#float PM17(Time, ROW) ;
PM17 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM17','f4',('Time','ROW'))
PM17[:] = TotlPoint_w_extra_00to12Z_sundy_PM17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM17'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM17'], varattr);
        setattr(PM17, varattr, varattrVal)
        
#float PM18(Time, ROW) ;
PM18 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM18','f4',('Time','ROW'))
PM18[:] = TotlPoint_w_extra_00to12Z_sundy_PM18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM18'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM18'], varattr);
        setattr(PM18, varattr, varattrVal)
        
#float PM19(Time, ROW) ;
PM19 = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('PM19','f4',('Time','ROW'))
PM19[:] = TotlPoint_w_extra_00to12Z_sundy_PM19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_00to12Z_sundy_file.variables['PM19'], varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file.variables['PM19'], varattr);
        setattr(PM19, varattr, varattrVal)

#char Times(Time) ;
Times = TotlPoint_w_extra_00to12Z_sundy_file.createVariable('Times','S1',('Time'))
Times[:] = TotlPoint_00to12Z_sundy_Times

#copy global attributes from TotlPoint_00to12Z_sundy_file
for varattr in TotlPoint_00to12Z_sundy_file.ncattrs():
    if hasattr(TotlPoint_00to12Z_sundy_file, varattr):
        varattrVal = getattr(TotlPoint_00to12Z_sundy_file, varattr);
        setattr(TotlPoint_w_extra_00to12Z_sundy_file, varattr, varattrVal)

TotlPoint_w_extra_00to12Z_sundy_file.close()


# In[21]:


#append extra points to original point file that is input to wrfchemi assembly program

###################################################################################################
#sundy, 12to24Z

###################################################################################################
#read original variables
TotlPoint_12to24Z_sundy_fn = base_dir+'/sundy/PtOG_12to24Z.nc'
TotlPoint_12to24Z_sundy_file = Dataset(TotlPoint_12to24Z_sundy_fn,mode='r',open=True)
TotlPoint_12to24Z_sundy_ITYPE = TotlPoint_12to24Z_sundy_file.variables['ITYPE'][:]
TotlPoint_12to24Z_sundy_STKht = TotlPoint_12to24Z_sundy_file.variables['STKht'][:]
TotlPoint_12to24Z_sundy_STKdiam = TotlPoint_12to24Z_sundy_file.variables['STKdiam'][:]
TotlPoint_12to24Z_sundy_STKtemp = TotlPoint_12to24Z_sundy_file.variables['STKtemp'][:]
TotlPoint_12to24Z_sundy_STKve = TotlPoint_12to24Z_sundy_file.variables['STKve'][:]
TotlPoint_12to24Z_sundy_STKflw = TotlPoint_12to24Z_sundy_file.variables['STKflw'][:]
TotlPoint_12to24Z_sundy_FUGht = TotlPoint_12to24Z_sundy_file.variables['FUGht'][:]
TotlPoint_12to24Z_sundy_XLONG = TotlPoint_12to24Z_sundy_file.variables['XLONG'][:]
TotlPoint_12to24Z_sundy_XLAT = TotlPoint_12to24Z_sundy_file.variables['XLAT'][:]
TotlPoint_12to24Z_sundy_CO2 = TotlPoint_12to24Z_sundy_file.variables['CO2'][:][:] 
TotlPoint_12to24Z_sundy_CO = TotlPoint_12to24Z_sundy_file.variables['CO'][:][:] 
TotlPoint_12to24Z_sundy_NH3 = TotlPoint_12to24Z_sundy_file.variables['NH3'][:][:] 
TotlPoint_12to24Z_sundy_NOX = TotlPoint_12to24Z_sundy_file.variables['NOX'][:][:] 
TotlPoint_12to24Z_sundy_PM10_PRI = TotlPoint_12to24Z_sundy_file.variables['PM10-PRI'][:][:] 
TotlPoint_12to24Z_sundy_PM25_PRI = TotlPoint_12to24Z_sundy_file.variables['PM25-PRI'][:][:] 
TotlPoint_12to24Z_sundy_SO2 = TotlPoint_12to24Z_sundy_file.variables['SO2'][:][:] 
TotlPoint_12to24Z_sundy_VOC = TotlPoint_12to24Z_sundy_file.variables['VOC'][:][:] 
TotlPoint_12to24Z_sundy_HC01 = TotlPoint_12to24Z_sundy_file.variables['HC01'][:][:] 
TotlPoint_12to24Z_sundy_HC02 = TotlPoint_12to24Z_sundy_file.variables['HC02'][:][:] 
TotlPoint_12to24Z_sundy_HC03 = TotlPoint_12to24Z_sundy_file.variables['HC03'][:][:] 
TotlPoint_12to24Z_sundy_HC04 = TotlPoint_12to24Z_sundy_file.variables['HC04'][:][:] 
TotlPoint_12to24Z_sundy_HC05 = TotlPoint_12to24Z_sundy_file.variables['HC05'][:][:] 
TotlPoint_12to24Z_sundy_HC06 = TotlPoint_12to24Z_sundy_file.variables['HC06'][:][:] 
TotlPoint_12to24Z_sundy_HC07 = TotlPoint_12to24Z_sundy_file.variables['HC07'][:][:] 
TotlPoint_12to24Z_sundy_HC08 = TotlPoint_12to24Z_sundy_file.variables['HC08'][:][:] 
TotlPoint_12to24Z_sundy_HC09 = TotlPoint_12to24Z_sundy_file.variables['HC09'][:][:] 
TotlPoint_12to24Z_sundy_HC10 = TotlPoint_12to24Z_sundy_file.variables['HC10'][:][:] 
TotlPoint_12to24Z_sundy_HC11 = TotlPoint_12to24Z_sundy_file.variables['HC11'][:][:] 
TotlPoint_12to24Z_sundy_HC12 = TotlPoint_12to24Z_sundy_file.variables['HC12'][:][:] 
TotlPoint_12to24Z_sundy_HC13 = TotlPoint_12to24Z_sundy_file.variables['HC13'][:][:] 
TotlPoint_12to24Z_sundy_HC14 = TotlPoint_12to24Z_sundy_file.variables['HC14'][:][:] 
TotlPoint_12to24Z_sundy_HC15 = TotlPoint_12to24Z_sundy_file.variables['HC15'][:][:] 
TotlPoint_12to24Z_sundy_HC16 = TotlPoint_12to24Z_sundy_file.variables['HC16'][:][:] 
TotlPoint_12to24Z_sundy_HC17 = TotlPoint_12to24Z_sundy_file.variables['HC17'][:][:] 
TotlPoint_12to24Z_sundy_HC18 = TotlPoint_12to24Z_sundy_file.variables['HC18'][:][:] 
TotlPoint_12to24Z_sundy_HC19 = TotlPoint_12to24Z_sundy_file.variables['HC19'][:][:] 
TotlPoint_12to24Z_sundy_HC20 = TotlPoint_12to24Z_sundy_file.variables['HC20'][:][:] 
TotlPoint_12to24Z_sundy_HC21 = TotlPoint_12to24Z_sundy_file.variables['HC21'][:][:] 
TotlPoint_12to24Z_sundy_HC22 = TotlPoint_12to24Z_sundy_file.variables['HC22'][:][:] 
TotlPoint_12to24Z_sundy_HC23 = TotlPoint_12to24Z_sundy_file.variables['HC23'][:][:] 
TotlPoint_12to24Z_sundy_HC24 = TotlPoint_12to24Z_sundy_file.variables['HC24'][:][:] 
TotlPoint_12to24Z_sundy_HC25 = TotlPoint_12to24Z_sundy_file.variables['HC25'][:][:] 
TotlPoint_12to24Z_sundy_HC26 = TotlPoint_12to24Z_sundy_file.variables['HC26'][:][:] 
TotlPoint_12to24Z_sundy_HC27 = TotlPoint_12to24Z_sundy_file.variables['HC27'][:][:] 
TotlPoint_12to24Z_sundy_HC28 = TotlPoint_12to24Z_sundy_file.variables['HC28'][:][:] 
TotlPoint_12to24Z_sundy_HC29 = TotlPoint_12to24Z_sundy_file.variables['HC29'][:][:] 
TotlPoint_12to24Z_sundy_HC30 = TotlPoint_12to24Z_sundy_file.variables['HC30'][:][:] 
TotlPoint_12to24Z_sundy_HC31 = TotlPoint_12to24Z_sundy_file.variables['HC31'][:][:] 
TotlPoint_12to24Z_sundy_HC32 = TotlPoint_12to24Z_sundy_file.variables['HC32'][:][:] 
TotlPoint_12to24Z_sundy_HC33 = TotlPoint_12to24Z_sundy_file.variables['HC33'][:][:] 
TotlPoint_12to24Z_sundy_HC34 = TotlPoint_12to24Z_sundy_file.variables['HC34'][:][:] 
TotlPoint_12to24Z_sundy_HC35 = TotlPoint_12to24Z_sundy_file.variables['HC35'][:][:] 
TotlPoint_12to24Z_sundy_HC36 = TotlPoint_12to24Z_sundy_file.variables['HC36'][:][:] 
TotlPoint_12to24Z_sundy_HC37 = TotlPoint_12to24Z_sundy_file.variables['HC37'][:][:] 
TotlPoint_12to24Z_sundy_HC38 = TotlPoint_12to24Z_sundy_file.variables['HC38'][:][:] 
TotlPoint_12to24Z_sundy_HC39 = TotlPoint_12to24Z_sundy_file.variables['HC39'][:][:] 
TotlPoint_12to24Z_sundy_HC40 = TotlPoint_12to24Z_sundy_file.variables['HC40'][:][:] 
TotlPoint_12to24Z_sundy_HC41 = TotlPoint_12to24Z_sundy_file.variables['HC41'][:][:] 
TotlPoint_12to24Z_sundy_HC42 = TotlPoint_12to24Z_sundy_file.variables['HC42'][:][:] 
TotlPoint_12to24Z_sundy_HC43 = TotlPoint_12to24Z_sundy_file.variables['HC43'][:][:] 
TotlPoint_12to24Z_sundy_HC44 = TotlPoint_12to24Z_sundy_file.variables['HC44'][:][:] 
TotlPoint_12to24Z_sundy_HC45 = TotlPoint_12to24Z_sundy_file.variables['HC45'][:][:] 
TotlPoint_12to24Z_sundy_HC46 = TotlPoint_12to24Z_sundy_file.variables['HC46'][:][:] 
TotlPoint_12to24Z_sundy_HC47 = TotlPoint_12to24Z_sundy_file.variables['HC47'][:][:] 
TotlPoint_12to24Z_sundy_HC48 = TotlPoint_12to24Z_sundy_file.variables['HC48'][:][:] 
TotlPoint_12to24Z_sundy_HC49 = TotlPoint_12to24Z_sundy_file.variables['HC49'][:][:] 
TotlPoint_12to24Z_sundy_HC50 = TotlPoint_12to24Z_sundy_file.variables['HC50'][:][:] 
TotlPoint_12to24Z_sundy_HC51 = TotlPoint_12to24Z_sundy_file.variables['HC51'][:][:] 
TotlPoint_12to24Z_sundy_HC52 = TotlPoint_12to24Z_sundy_file.variables['HC52'][:][:] 
TotlPoint_12to24Z_sundy_HC53 = TotlPoint_12to24Z_sundy_file.variables['HC53'][:][:] 
TotlPoint_12to24Z_sundy_HC54 = TotlPoint_12to24Z_sundy_file.variables['HC54'][:][:] 
TotlPoint_12to24Z_sundy_HC55 = TotlPoint_12to24Z_sundy_file.variables['HC55'][:][:] 
TotlPoint_12to24Z_sundy_HC56 = TotlPoint_12to24Z_sundy_file.variables['HC56'][:][:] 
TotlPoint_12to24Z_sundy_HC57 = TotlPoint_12to24Z_sundy_file.variables['HC57'][:][:] 
TotlPoint_12to24Z_sundy_HC58 = TotlPoint_12to24Z_sundy_file.variables['HC58'][:][:] 
TotlPoint_12to24Z_sundy_HC59 = TotlPoint_12to24Z_sundy_file.variables['HC59'][:][:] 
TotlPoint_12to24Z_sundy_HC60 = TotlPoint_12to24Z_sundy_file.variables['HC60'][:][:] 
TotlPoint_12to24Z_sundy_HC61 = TotlPoint_12to24Z_sundy_file.variables['HC61'][:][:] 
TotlPoint_12to24Z_sundy_HC62 = TotlPoint_12to24Z_sundy_file.variables['HC62'][:][:] 
TotlPoint_12to24Z_sundy_HC63 = TotlPoint_12to24Z_sundy_file.variables['HC63'][:][:] 
TotlPoint_12to24Z_sundy_HC64 = TotlPoint_12to24Z_sundy_file.variables['HC64'][:][:] 
TotlPoint_12to24Z_sundy_HC65 = TotlPoint_12to24Z_sundy_file.variables['HC65'][:][:] 
TotlPoint_12to24Z_sundy_HC66 = TotlPoint_12to24Z_sundy_file.variables['HC66'][:][:] 
TotlPoint_12to24Z_sundy_HC67 = TotlPoint_12to24Z_sundy_file.variables['HC67'][:][:] 
TotlPoint_12to24Z_sundy_HC68 = TotlPoint_12to24Z_sundy_file.variables['HC68'][:][:] 
TotlPoint_12to24Z_sundy_PM01 = TotlPoint_12to24Z_sundy_file.variables['PM01'][:][:] 
TotlPoint_12to24Z_sundy_PM02 = TotlPoint_12to24Z_sundy_file.variables['PM02'][:][:] 
TotlPoint_12to24Z_sundy_PM03 = TotlPoint_12to24Z_sundy_file.variables['PM03'][:][:] 
TotlPoint_12to24Z_sundy_PM04 = TotlPoint_12to24Z_sundy_file.variables['PM04'][:][:] 
TotlPoint_12to24Z_sundy_PM05 = TotlPoint_12to24Z_sundy_file.variables['PM05'][:][:] 
TotlPoint_12to24Z_sundy_PM06 = TotlPoint_12to24Z_sundy_file.variables['PM06'][:][:] 
TotlPoint_12to24Z_sundy_PM07 = TotlPoint_12to24Z_sundy_file.variables['PM07'][:][:] 
TotlPoint_12to24Z_sundy_PM08 = TotlPoint_12to24Z_sundy_file.variables['PM08'][:][:] 
TotlPoint_12to24Z_sundy_PM09 = TotlPoint_12to24Z_sundy_file.variables['PM09'][:][:] 
TotlPoint_12to24Z_sundy_PM10 = TotlPoint_12to24Z_sundy_file.variables['PM10'][:][:] 
TotlPoint_12to24Z_sundy_PM11 = TotlPoint_12to24Z_sundy_file.variables['PM11'][:][:] 
TotlPoint_12to24Z_sundy_PM12 = TotlPoint_12to24Z_sundy_file.variables['PM12'][:][:] 
TotlPoint_12to24Z_sundy_PM13 = TotlPoint_12to24Z_sundy_file.variables['PM13'][:][:] 
TotlPoint_12to24Z_sundy_PM14 = TotlPoint_12to24Z_sundy_file.variables['PM14'][:][:] 
TotlPoint_12to24Z_sundy_PM15 = TotlPoint_12to24Z_sundy_file.variables['PM15'][:][:] 
TotlPoint_12to24Z_sundy_PM16 = TotlPoint_12to24Z_sundy_file.variables['PM16'][:][:] 
TotlPoint_12to24Z_sundy_PM17 = TotlPoint_12to24Z_sundy_file.variables['PM17'][:][:] 
TotlPoint_12to24Z_sundy_PM18 = TotlPoint_12to24Z_sundy_file.variables['PM18'][:][:] 
TotlPoint_12to24Z_sundy_PM19 = TotlPoint_12to24Z_sundy_file.variables['PM19'][:][:] 
TotlPoint_12to24Z_sundy_Times = TotlPoint_12to24Z_sundy_file.variables['Times'][:]

###################################################################################################
#get total ROW
nROW_org, = TotlPoint_12to24Z_sundy_ITYPE.shape
nROW_extra_OG = len(LON_ng_proc)
nROW_extra = nROW_extra_OG
nROW = nROW_org + nROW_extra
print("nROW_org",nROW_org)
print("nROW_extra",nROW_extra)
print("nROW",nROW)

###################################################################################################
#Organize extra_data
extra_ITYPE_OG = np.array(ERPTYPE_ng_proc)
extra_ITYPE = extra_ITYPE_OG

extra_STKht_OG = np.array(STKHGT_ng_proc)
extra_STKht = extra_STKht_OG

extra_STKdiam_OG = np.array(STKDIAM_ng_proc)
extra_STKdiam = extra_STKdiam_OG

extra_STKtemp_OG = np.array(STKTEMP_ng_proc)
extra_STKtemp = extra_STKtemp_OG

extra_STKve_OG = np.array(STKVEL_ng_proc)
extra_STKve = extra_STKve_OG

extra_STKflw_OG = np.array(STKFLOW_ng_proc)
extra_STKflw = extra_STKflw_OG

extra_FUGht = np.empty(nROW_extra) #FUGht set as empty

extra_XLONG_OG = np.array(LON_ng_proc)
extra_XLONG = extra_XLONG_OG

extra_XLAT_OG = np.array(LAT_ng_proc)
extra_XLAT = extra_XLAT_OG

extra_STATE_OG = STATE_ng_proc

###################################################################################################
#CO2

##################################################################################
extra_CO2_FCPE_ng_proc = HRall_CO2_FCPE_MetricTon_2022mm_ng_proc_sundy[12:24,:]
extra_CO2_OG = extra_CO2_FCPE_ng_proc

##################################################################################
extra_CO2 =extra_CO2_OG

###################################################################################################
#CH4 from OG can use GHGRP numbers

##################################################################################
extra_CH4_FCPE_ng_proc = HRall_CH4_FCPE_MetricTon_2022mm_ng_proc_sundy[12:24,:]
extra_CH4_OG = extra_CH4_FCPE_ng_proc

###################################################################################################

species_vector = ['CO','NH3','NOX','PM10-PRI','PM25-PRI','SO2','VOC',
                  'HC01','HC02','HC03','HC04','HC05','HC06','HC07','HC08','HC09','HC10',
                  'HC11','HC12','HC13','HC14','HC15','HC16','HC17','HC18','HC19','HC20',
                  'HC21','HC22','HC23','HC24','HC25','HC26','HC27','HC28','HC29','HC30',
                  'HC31','HC32','HC33','HC34','HC35','HC36','HC37','HC38','HC39','HC40',
                  'HC41','HC42','HC43','HC44','HC45','HC46','HC47','HC48','HC49','HC50',
                  'PM01','PM02','PM03','PM04','PM05','PM06','PM07','PM08','PM09','PM10',
                  'PM11','PM12','PM13','PM14','PM15','PM16','PM17','PM18','PM19']

states_vector = ['Alabama','Arizona','Arkansas','California','Colorado','Connecticut',
                 'Delaware','District of Columbia','Florida','Georgia','Idaho','Illinois','Indiana','Iowa',
                 'Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts',
                 'Michigan','Minnesota','Mississippi','Missouri','Montana','Nebraska',
                 'Nevada','New Hampshire','New Jersey','New Mexico','New York',
                 'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania',
                 'Rhode Island','South Carolina','South Dakota','Tennessee','Texas','Utah',
                 'Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

states_abb_vector = ['AL', 'AZ', 'AR', 'CA', 'CO', 'CT', 
                     'DE', 'DC', 'FL', 'GA', 'ID', 'IL', 'IN', 'IA', 
                     'KS', 'KY', 'LA', 'ME', 'MD', 'MA', 
                     'MI', 'MN', 'MS', 'MO', 'MT', 'NE', 
                     'NV', 'NH', 'NJ', 'NM', 'NY', 
                     'NC', 'ND', 'OH', 'OK', 'OR', 'PA', 
                     'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 
                     'VT', 'VA', 'WA', 'WV', 'WI', 'WY']

###################################################################################################
#AQ: using state-level emission ratios to CO2

###################################################################################################
#grab emission ratios from the summary ratio arrays
###################################################################################################

#based on OG process type, species, and state location
extra_X_OG = np.empty([12,nROW_extra_OG,len(species_vector)])
print("extra_X_OG.shape", extra_X_OG.shape)

for pt in range(0,nROW_extra_OG):
    proc_index = 0
    #print("proc_index",proc_index)

    state_abb_cur = extra_STATE_OG[pt]
    
    if state_abb_cur in states_abb_vector:
        state_index = states_abb_vector.index(state_abb_cur)
        #print("state_index",state_index)
    
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * proc_spec_state_emisXdCO2_OG[proc_index,spec_index,state_index]

    else:
        for spec_cur in species_vector:
            spec_index = species_vector.index(spec_cur)
            #print("spec_index",spec_index)
            extra_X_OG[:,pt,spec_index] = extra_CO2_OG[:,pt] * statistics.mean(proc_spec_state_emisXdCO2_OG[proc_index,spec_index,:])

extra_X_OG_dict = {}
for spec_cur in species_vector:
    #for CH4 use GHGRP numbers
    if spec_cur == 'HC01':
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_CH4_OG
    else:
        spec_index = species_vector.index(spec_cur)
        #print("spec_index",spec_index)
        extra_X_OG_dict["extra_{0}_OG".format(spec_cur)] = extra_X_OG[:,:,spec_index]

###################################################################################################
#stack OG AQ species
extra_X_dict = {}
for spec_cur in species_vector:
    extra_Xi_OG = extra_X_OG_dict["extra_{0}_OG".format(spec_cur)]
    extra_X = extra_Xi_OG
    extra_X_dict["extra_{0}".format(spec_cur)] = extra_X

###################################################################################################
extra_other_spec = np.zeros([12,nROW_extra])

###################################################################################################
#append extra points to original data
TotlPoint_w_extra_12to24Z_sundy_ITYPE = np.concatenate((TotlPoint_12to24Z_sundy_ITYPE,extra_ITYPE),axis=0)
TotlPoint_w_extra_12to24Z_sundy_STKht = np.concatenate((TotlPoint_12to24Z_sundy_STKht,extra_STKht),axis=0)
TotlPoint_w_extra_12to24Z_sundy_STKdiam = np.concatenate((TotlPoint_12to24Z_sundy_STKdiam,extra_STKdiam),axis=0)
TotlPoint_w_extra_12to24Z_sundy_STKtemp = np.concatenate((TotlPoint_12to24Z_sundy_STKtemp,extra_STKtemp),axis=0)
TotlPoint_w_extra_12to24Z_sundy_STKve = np.concatenate((TotlPoint_12to24Z_sundy_STKve,extra_STKve),axis=0)
TotlPoint_w_extra_12to24Z_sundy_STKflw = np.concatenate((TotlPoint_12to24Z_sundy_STKflw,extra_STKflw),axis=0)
TotlPoint_w_extra_12to24Z_sundy_FUGht = np.concatenate((TotlPoint_12to24Z_sundy_FUGht,extra_FUGht),axis=0)
TotlPoint_w_extra_12to24Z_sundy_XLONG = np.concatenate((TotlPoint_12to24Z_sundy_XLONG,extra_XLONG),axis=0)
TotlPoint_w_extra_12to24Z_sundy_XLAT = np.concatenate((TotlPoint_12to24Z_sundy_XLAT,extra_XLAT),axis=0)
TotlPoint_w_extra_12to24Z_sundy_CO2 = np.concatenate((TotlPoint_12to24Z_sundy_CO2,extra_CO2),axis=1)
TotlPoint_w_extra_12to24Z_sundy_CO = np.concatenate((TotlPoint_12to24Z_sundy_CO,extra_X_dict["extra_CO"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_NH3 = np.concatenate((TotlPoint_12to24Z_sundy_NH3,extra_X_dict["extra_NH3"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_NOX = np.concatenate((TotlPoint_12to24Z_sundy_NOX,extra_X_dict["extra_NOX"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM10_PRI = np.concatenate((TotlPoint_12to24Z_sundy_PM10_PRI,extra_X_dict["extra_PM10-PRI"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM25_PRI = np.concatenate((TotlPoint_12to24Z_sundy_PM25_PRI,extra_X_dict["extra_PM25-PRI"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_SO2 = np.concatenate((TotlPoint_12to24Z_sundy_SO2,extra_X_dict["extra_SO2"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_VOC = np.concatenate((TotlPoint_12to24Z_sundy_VOC,extra_X_dict["extra_VOC"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC01 = np.concatenate((TotlPoint_12to24Z_sundy_HC01,extra_X_dict["extra_HC01"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC02 = np.concatenate((TotlPoint_12to24Z_sundy_HC02,extra_X_dict["extra_HC02"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC03 = np.concatenate((TotlPoint_12to24Z_sundy_HC03,extra_X_dict["extra_HC03"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC04 = np.concatenate((TotlPoint_12to24Z_sundy_HC04,extra_X_dict["extra_HC04"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC05 = np.concatenate((TotlPoint_12to24Z_sundy_HC05,extra_X_dict["extra_HC05"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC06 = np.concatenate((TotlPoint_12to24Z_sundy_HC06,extra_X_dict["extra_HC06"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC07 = np.concatenate((TotlPoint_12to24Z_sundy_HC07,extra_X_dict["extra_HC07"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC08 = np.concatenate((TotlPoint_12to24Z_sundy_HC08,extra_X_dict["extra_HC08"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC09 = np.concatenate((TotlPoint_12to24Z_sundy_HC09,extra_X_dict["extra_HC09"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC10 = np.concatenate((TotlPoint_12to24Z_sundy_HC10,extra_X_dict["extra_HC10"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC11 = np.concatenate((TotlPoint_12to24Z_sundy_HC11,extra_X_dict["extra_HC11"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC12 = np.concatenate((TotlPoint_12to24Z_sundy_HC12,extra_X_dict["extra_HC12"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC13 = np.concatenate((TotlPoint_12to24Z_sundy_HC13,extra_X_dict["extra_HC13"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC14 = np.concatenate((TotlPoint_12to24Z_sundy_HC14,extra_X_dict["extra_HC14"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC15 = np.concatenate((TotlPoint_12to24Z_sundy_HC15,extra_X_dict["extra_HC15"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC16 = np.concatenate((TotlPoint_12to24Z_sundy_HC16,extra_X_dict["extra_HC16"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC17 = np.concatenate((TotlPoint_12to24Z_sundy_HC17,extra_X_dict["extra_HC17"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC18 = np.concatenate((TotlPoint_12to24Z_sundy_HC18,extra_X_dict["extra_HC18"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC19 = np.concatenate((TotlPoint_12to24Z_sundy_HC19,extra_X_dict["extra_HC19"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC20 = np.concatenate((TotlPoint_12to24Z_sundy_HC20,extra_X_dict["extra_HC20"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC21 = np.concatenate((TotlPoint_12to24Z_sundy_HC21,extra_X_dict["extra_HC21"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC22 = np.concatenate((TotlPoint_12to24Z_sundy_HC22,extra_X_dict["extra_HC22"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC23 = np.concatenate((TotlPoint_12to24Z_sundy_HC23,extra_X_dict["extra_HC23"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC24 = np.concatenate((TotlPoint_12to24Z_sundy_HC24,extra_X_dict["extra_HC24"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC25 = np.concatenate((TotlPoint_12to24Z_sundy_HC25,extra_X_dict["extra_HC25"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC26 = np.concatenate((TotlPoint_12to24Z_sundy_HC26,extra_X_dict["extra_HC26"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC27 = np.concatenate((TotlPoint_12to24Z_sundy_HC27,extra_X_dict["extra_HC27"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC28 = np.concatenate((TotlPoint_12to24Z_sundy_HC28,extra_X_dict["extra_HC28"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC29 = np.concatenate((TotlPoint_12to24Z_sundy_HC29,extra_X_dict["extra_HC29"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC30 = np.concatenate((TotlPoint_12to24Z_sundy_HC30,extra_X_dict["extra_HC30"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC31 = np.concatenate((TotlPoint_12to24Z_sundy_HC31,extra_X_dict["extra_HC31"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC32 = np.concatenate((TotlPoint_12to24Z_sundy_HC32,extra_X_dict["extra_HC32"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC33 = np.concatenate((TotlPoint_12to24Z_sundy_HC33,extra_X_dict["extra_HC33"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC34 = np.concatenate((TotlPoint_12to24Z_sundy_HC34,extra_X_dict["extra_HC34"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC35 = np.concatenate((TotlPoint_12to24Z_sundy_HC35,extra_X_dict["extra_HC35"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC36 = np.concatenate((TotlPoint_12to24Z_sundy_HC36,extra_X_dict["extra_HC36"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC37 = np.concatenate((TotlPoint_12to24Z_sundy_HC37,extra_X_dict["extra_HC37"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC38 = np.concatenate((TotlPoint_12to24Z_sundy_HC38,extra_X_dict["extra_HC38"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC39 = np.concatenate((TotlPoint_12to24Z_sundy_HC39,extra_X_dict["extra_HC39"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC40 = np.concatenate((TotlPoint_12to24Z_sundy_HC40,extra_X_dict["extra_HC40"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC41 = np.concatenate((TotlPoint_12to24Z_sundy_HC41,extra_X_dict["extra_HC41"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC42 = np.concatenate((TotlPoint_12to24Z_sundy_HC42,extra_X_dict["extra_HC42"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC43 = np.concatenate((TotlPoint_12to24Z_sundy_HC43,extra_X_dict["extra_HC43"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC44 = np.concatenate((TotlPoint_12to24Z_sundy_HC44,extra_X_dict["extra_HC44"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC45 = np.concatenate((TotlPoint_12to24Z_sundy_HC45,extra_X_dict["extra_HC45"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC46 = np.concatenate((TotlPoint_12to24Z_sundy_HC46,extra_X_dict["extra_HC46"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC47 = np.concatenate((TotlPoint_12to24Z_sundy_HC47,extra_X_dict["extra_HC47"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC48 = np.concatenate((TotlPoint_12to24Z_sundy_HC48,extra_X_dict["extra_HC48"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC49 = np.concatenate((TotlPoint_12to24Z_sundy_HC49,extra_X_dict["extra_HC49"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC50 = np.concatenate((TotlPoint_12to24Z_sundy_HC50,extra_X_dict["extra_HC50"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC51 = np.concatenate((TotlPoint_12to24Z_sundy_HC51,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC52 = np.concatenate((TotlPoint_12to24Z_sundy_HC52,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC53 = np.concatenate((TotlPoint_12to24Z_sundy_HC53,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC54 = np.concatenate((TotlPoint_12to24Z_sundy_HC54,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC55 = np.concatenate((TotlPoint_12to24Z_sundy_HC55,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC56 = np.concatenate((TotlPoint_12to24Z_sundy_HC56,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC57 = np.concatenate((TotlPoint_12to24Z_sundy_HC57,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC58 = np.concatenate((TotlPoint_12to24Z_sundy_HC58,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC59 = np.concatenate((TotlPoint_12to24Z_sundy_HC59,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC60 = np.concatenate((TotlPoint_12to24Z_sundy_HC60,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC61 = np.concatenate((TotlPoint_12to24Z_sundy_HC61,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC62 = np.concatenate((TotlPoint_12to24Z_sundy_HC62,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC63 = np.concatenate((TotlPoint_12to24Z_sundy_HC63,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC64 = np.concatenate((TotlPoint_12to24Z_sundy_HC64,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC65 = np.concatenate((TotlPoint_12to24Z_sundy_HC65,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC66 = np.concatenate((TotlPoint_12to24Z_sundy_HC66,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC67 = np.concatenate((TotlPoint_12to24Z_sundy_HC67,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_HC68 = np.concatenate((TotlPoint_12to24Z_sundy_HC68,extra_other_spec),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM01 = np.concatenate((TotlPoint_12to24Z_sundy_PM01,extra_X_dict["extra_PM01"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM02 = np.concatenate((TotlPoint_12to24Z_sundy_PM02,extra_X_dict["extra_PM02"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM03 = np.concatenate((TotlPoint_12to24Z_sundy_PM03,extra_X_dict["extra_PM03"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM04 = np.concatenate((TotlPoint_12to24Z_sundy_PM04,extra_X_dict["extra_PM04"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM05 = np.concatenate((TotlPoint_12to24Z_sundy_PM05,extra_X_dict["extra_PM05"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM06 = np.concatenate((TotlPoint_12to24Z_sundy_PM06,extra_X_dict["extra_PM06"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM07 = np.concatenate((TotlPoint_12to24Z_sundy_PM07,extra_X_dict["extra_PM07"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM08 = np.concatenate((TotlPoint_12to24Z_sundy_PM08,extra_X_dict["extra_PM08"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM09 = np.concatenate((TotlPoint_12to24Z_sundy_PM09,extra_X_dict["extra_PM09"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM10 = np.concatenate((TotlPoint_12to24Z_sundy_PM10,extra_X_dict["extra_PM10"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM11 = np.concatenate((TotlPoint_12to24Z_sundy_PM11,extra_X_dict["extra_PM11"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM12 = np.concatenate((TotlPoint_12to24Z_sundy_PM12,extra_X_dict["extra_PM12"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM13 = np.concatenate((TotlPoint_12to24Z_sundy_PM13,extra_X_dict["extra_PM13"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM14 = np.concatenate((TotlPoint_12to24Z_sundy_PM14,extra_X_dict["extra_PM14"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM15 = np.concatenate((TotlPoint_12to24Z_sundy_PM15,extra_X_dict["extra_PM15"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM16 = np.concatenate((TotlPoint_12to24Z_sundy_PM16,extra_X_dict["extra_PM16"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM17 = np.concatenate((TotlPoint_12to24Z_sundy_PM17,extra_X_dict["extra_PM17"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM18 = np.concatenate((TotlPoint_12to24Z_sundy_PM18,extra_X_dict["extra_PM18"]),axis=1)
TotlPoint_w_extra_12to24Z_sundy_PM19 = np.concatenate((TotlPoint_12to24Z_sundy_PM19,extra_X_dict["extra_PM19"]),axis=1)

###################################################################################################
#write total points with extra points appended
TotlPoint_w_extra_12to24Z_sundy_fn = append_dir+'/sundy/PtOG_12to24Z.nc'
TotlPoint_w_extra_12to24Z_sundy_file = Dataset(TotlPoint_w_extra_12to24Z_sundy_fn,mode='w',format='NETCDF3_64BIT')

#Creat dimensions
TotlPoint_w_extra_12to24Z_sundy_file.createDimension("ROW", nROW)
TotlPoint_w_extra_12to24Z_sundy_file.createDimension("Time", 12)
TotlPoint_w_extra_12to24Z_sundy_file.sync()

#Create variables
#float ITYPE(ROW) ;
ITYPE = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('ITYPE','f4',('ROW'),fill_value = float(0))
ITYPE[:] = TotlPoint_w_extra_12to24Z_sundy_ITYPE
varattrs=["FieldType","MemoryOrder","description","units","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['ITYPE'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['ITYPE'], varattr);
        setattr(ITYPE, varattr, varattrVal)

#float STKht(ROW) ;
STKht = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('STKht','f4',('ROW'),fill_value = float(0))
STKht[:] = TotlPoint_w_extra_12to24Z_sundy_STKht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['STKht'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['STKht'], varattr);
        setattr(STKht, varattr, varattrVal)
        
#float STKdiam(ROW) ;
STKdiam = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('STKdiam','f4',('ROW'),fill_value = float(0))
STKdiam[:] = TotlPoint_w_extra_12to24Z_sundy_STKdiam
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['STKdiam'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['STKdiam'], varattr);
        setattr(STKdiam, varattr, varattrVal)

#float STKtemp(ROW) ;
STKtemp = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('STKtemp','f4',('ROW'),fill_value = float(0))
STKtemp[:] = TotlPoint_w_extra_12to24Z_sundy_STKtemp
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['STKtemp'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['STKtemp'], varattr);
        setattr(STKtemp, varattr, varattrVal)
        
#float STKve(ROW) ;
STKve = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('STKve','f4',('ROW'),fill_value = float(0))
STKve[:] = TotlPoint_w_extra_12to24Z_sundy_STKve
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['STKve'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['STKve'], varattr);
        setattr(STKve, varattr, varattrVal)
        
#float STKflw(ROW) ;
STKflw = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('STKflw','f4',('ROW'),fill_value = float(0))
STKflw[:] = TotlPoint_w_extra_12to24Z_sundy_STKflw
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['STKflw'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['STKflw'], varattr);
        setattr(STKflw, varattr, varattrVal)
        
#float FUGht(ROW) ;
FUGht = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('FUGht','f4',('ROW'),fill_value = float(0))
FUGht[:] = TotlPoint_w_extra_12to24Z_sundy_FUGht
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['FUGht'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['FUGht'], varattr);
        setattr(FUGht, varattr, varattrVal)
        
#float XLONG(ROW) ;
XLONG = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('XLONG','f4',('ROW'),fill_value = float(0))
XLONG[:] = TotlPoint_w_extra_12to24Z_sundy_XLONG
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['XLONG'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['XLONG'], varattr);
        setattr(XLONG, varattr, varattrVal)
        
#float XLAT(ROW) ;
XLAT = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('XLAT','f4',('ROW'),fill_value = float(0))
XLAT[:] = TotlPoint_w_extra_12to24Z_sundy_XLAT
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['XLAT'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['XLAT'], varattr);
        setattr(XLAT, varattr, varattrVal)

#float CO2(Time, ROW) ;
CO2 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('CO2','f4',('Time','ROW'))
CO2[:] = TotlPoint_w_extra_12to24Z_sundy_CO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['CO2'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['CO2'], varattr);
        setattr(CO2, varattr, varattrVal)
        
#float CO(Time, ROW) ;
CO = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('CO','f4',('Time','ROW'))
CO[:] = TotlPoint_w_extra_12to24Z_sundy_CO
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['CO'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['CO'], varattr);
        setattr(CO, varattr, varattrVal)
        
#float NH3(Time, ROW) ;
NH3 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('NH3','f4',('Time','ROW'))
NH3[:] = TotlPoint_w_extra_12to24Z_sundy_NH3
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['NH3'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['NH3'], varattr);
        setattr(NH3, varattr, varattrVal)
        
#float NOX(Time, ROW) ;
NOX = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('NOX','f4',('Time','ROW'))
NOX[:] = TotlPoint_w_extra_12to24Z_sundy_NOX
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['NOX'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['NOX'], varattr);
        setattr(NOX, varattr, varattrVal)
        
#float PM10-PRI(Time, ROW) ;
PM10_PRI = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM10-PRI','f4',('Time','ROW'))
PM10_PRI[:] = TotlPoint_w_extra_12to24Z_sundy_PM10_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM10-PRI'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM10-PRI'], varattr);
        setattr(PM10_PRI, varattr, varattrVal)
        
#float PM25-PRI(Time, ROW) ;
PM25_PRI = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM25-PRI','f4',('Time','ROW'))
PM25_PRI[:] = TotlPoint_w_extra_12to24Z_sundy_PM25_PRI
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM25-PRI'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM25-PRI'], varattr);
        setattr(PM25_PRI, varattr, varattrVal)
        
#float SO2(Time, ROW) ;
SO2 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('SO2','f4',('Time','ROW'))
SO2[:] = TotlPoint_w_extra_12to24Z_sundy_SO2
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['SO2'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['SO2'], varattr);
        setattr(SO2, varattr, varattrVal)
        
#float VOC(Time, ROW) ;
VOC = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('VOC','f4',('Time','ROW'))
VOC[:] = TotlPoint_w_extra_12to24Z_sundy_VOC
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['VOC'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['VOC'], varattr);
        setattr(VOC, varattr, varattrVal)
        
#float HC01(Time, ROW) ;
HC01 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC01','f4',('Time','ROW'))
HC01[:] = TotlPoint_w_extra_12to24Z_sundy_HC01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC01'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC01'], varattr);
        setattr(HC01, varattr, varattrVal)
        
#float HC02(Time, ROW) ;
HC02 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC02','f4',('Time','ROW'))
HC02[:] = TotlPoint_w_extra_12to24Z_sundy_HC02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC02'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC02'], varattr);
        setattr(HC02, varattr, varattrVal)
        
#float HC03(Time, ROW) ;
HC03 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC03','f4',('Time','ROW'))
HC03[:] = TotlPoint_w_extra_12to24Z_sundy_HC03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC03'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC03'], varattr);
        setattr(HC03, varattr, varattrVal)
        
#float HC04(Time, ROW) ;
HC04 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC04','f4',('Time','ROW'))
HC04[:] = TotlPoint_w_extra_12to24Z_sundy_HC04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC04'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC04'], varattr);
        setattr(HC04, varattr, varattrVal)
        
#float HC05(Time, ROW) ;
HC05 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC05','f4',('Time','ROW'))
HC05[:] = TotlPoint_w_extra_12to24Z_sundy_HC05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC05'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC05'], varattr);
        setattr(HC05, varattr, varattrVal)
        
#float HC06(Time, ROW) ;
HC06 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC06','f4',('Time','ROW'))
HC06[:] = TotlPoint_w_extra_12to24Z_sundy_HC06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC06'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC06'], varattr);
        setattr(HC06, varattr, varattrVal)
        
#float HC07(Time, ROW) ;
HC07 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC07','f4',('Time','ROW'))
HC07[:] = TotlPoint_w_extra_12to24Z_sundy_HC07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC07'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC07'], varattr);
        setattr(HC07, varattr, varattrVal)
        
#float HC08(Time, ROW) ;
HC08 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC08','f4',('Time','ROW'))
HC08[:] = TotlPoint_w_extra_12to24Z_sundy_HC08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC08'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC08'], varattr);
        setattr(HC08, varattr, varattrVal)
        
#float HC09(Time, ROW) ;
HC09 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC09','f4',('Time','ROW'))
HC09[:] = TotlPoint_w_extra_12to24Z_sundy_HC09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC09'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC09'], varattr);
        setattr(HC09, varattr, varattrVal)
        
#float HC10(Time, ROW) ;
HC10 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC10','f4',('Time','ROW'))
HC10[:] = TotlPoint_w_extra_12to24Z_sundy_HC10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC10'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC10'], varattr);
        setattr(HC10, varattr, varattrVal)

#float HC11(Time, ROW) ;
HC11 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC11','f4',('Time','ROW'))
HC11[:] = TotlPoint_w_extra_12to24Z_sundy_HC11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC11'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC11'], varattr);
        setattr(HC11, varattr, varattrVal)
        
#float HC12(Time, ROW) ;
HC12 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC12','f4',('Time','ROW'))
HC12[:] = TotlPoint_w_extra_12to24Z_sundy_HC12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC12'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC12'], varattr);
        setattr(HC12, varattr, varattrVal)
        
#float HC13(Time, ROW) ;
HC13 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC13','f4',('Time','ROW'))
HC13[:] = TotlPoint_w_extra_12to24Z_sundy_HC13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC13'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC13'], varattr);
        setattr(HC13, varattr, varattrVal)
        
#float HC14(Time, ROW) ;
HC14 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC14','f4',('Time','ROW'))
HC14[:] = TotlPoint_w_extra_12to24Z_sundy_HC14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC14'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC14'], varattr);
        setattr(HC14, varattr, varattrVal)
        
#float HC15(Time, ROW) ;
HC15 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC15','f4',('Time','ROW'))
HC15[:] = TotlPoint_w_extra_12to24Z_sundy_HC15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC15'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC15'], varattr);
        setattr(HC15, varattr, varattrVal)
        
#float HC16(Time, ROW) ;
HC16 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC16','f4',('Time','ROW'))
HC16[:] = TotlPoint_w_extra_12to24Z_sundy_HC16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC16'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC16'], varattr);
        setattr(HC16, varattr, varattrVal)
        
#float HC17(Time, ROW) ;
HC17 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC17','f4',('Time','ROW'))
HC17[:] = TotlPoint_w_extra_12to24Z_sundy_HC17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC17'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC17'], varattr);
        setattr(HC17, varattr, varattrVal)
        
#float HC18(Time, ROW) ;
HC18 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC18','f4',('Time','ROW'))
HC18[:] = TotlPoint_w_extra_12to24Z_sundy_HC18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC18'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC18'], varattr);
        setattr(HC18, varattr, varattrVal)
        
#float HC19(Time, ROW) ;
HC19 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC19','f4',('Time','ROW'))
HC19[:] = TotlPoint_w_extra_12to24Z_sundy_HC19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC19'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC19'], varattr);
        setattr(HC19, varattr, varattrVal)

#float HC20(Time, ROW) ;
HC20 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC20','f4',('Time','ROW'))
HC20[:] = TotlPoint_w_extra_12to24Z_sundy_HC20
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC20'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC20'], varattr);
        setattr(HC20, varattr, varattrVal)

#float HC21(Time, ROW) ;
HC21 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC21','f4',('Time','ROW'))
HC21[:] = TotlPoint_w_extra_12to24Z_sundy_HC21
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC21'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC21'], varattr);
        setattr(HC21, varattr, varattrVal)
        
#float HC22(Time, ROW) ;
HC22 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC22','f4',('Time','ROW'))
HC22[:] = TotlPoint_w_extra_12to24Z_sundy_HC22
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC22'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC22'], varattr);
        setattr(HC22, varattr, varattrVal)
        
#float HC23(Time, ROW) ;
HC23 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC23','f4',('Time','ROW'))
HC23[:] = TotlPoint_w_extra_12to24Z_sundy_HC23
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC23'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC23'], varattr);
        setattr(HC23, varattr, varattrVal)
        
#float HC24(Time, ROW) ;
HC24 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC24','f4',('Time','ROW'))
HC24[:] = TotlPoint_w_extra_12to24Z_sundy_HC24
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC24'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC24'], varattr);
        setattr(HC24, varattr, varattrVal)
        
#float HC25(Time, ROW) ;
HC25 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC25','f4',('Time','ROW'))
HC25[:] = TotlPoint_w_extra_12to24Z_sundy_HC25
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC25'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC25'], varattr);
        setattr(HC25, varattr, varattrVal)
        
#float HC26(Time, ROW) ;
HC26 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC26','f4',('Time','ROW'))
HC26[:] = TotlPoint_w_extra_12to24Z_sundy_HC26
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC26'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC26'], varattr);
        setattr(HC26, varattr, varattrVal)
        
#float HC27(Time, ROW) ;
HC27 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC27','f4',('Time','ROW'))
HC27[:] = TotlPoint_w_extra_12to24Z_sundy_HC27
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC27'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC27'], varattr);
        setattr(HC27, varattr, varattrVal)
        
#float HC28(Time, ROW) ;
HC28 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC28','f4',('Time','ROW'))
HC28[:] = TotlPoint_w_extra_12to24Z_sundy_HC28
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC28'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC28'], varattr);
        setattr(HC28, varattr, varattrVal)
        
#float HC29(Time, ROW) ;
HC29 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC29','f4',('Time','ROW'))
HC29[:] = TotlPoint_w_extra_12to24Z_sundy_HC29
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC29'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC29'], varattr);
        setattr(HC29, varattr, varattrVal)

#float HC30(Time, ROW) ;
HC30 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC30','f4',('Time','ROW'))
HC30[:] = TotlPoint_w_extra_12to24Z_sundy_HC30
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC30'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC30'], varattr);
        setattr(HC30, varattr, varattrVal)

#float HC31(Time, ROW) ;
HC31 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC31','f4',('Time','ROW'))
HC31[:] = TotlPoint_w_extra_12to24Z_sundy_HC31
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC31'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC31'], varattr);
        setattr(HC31, varattr, varattrVal)
        
#float HC32(Time, ROW) ;
HC32 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC32','f4',('Time','ROW'))
HC32[:] = TotlPoint_w_extra_12to24Z_sundy_HC32
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC32'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC32'], varattr);
        setattr(HC32, varattr, varattrVal)
        
#float HC33(Time, ROW) ;
HC33 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC33','f4',('Time','ROW'))
HC33[:] = TotlPoint_w_extra_12to24Z_sundy_HC33
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC33'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC33'], varattr);
        setattr(HC33, varattr, varattrVal)
        
#float HC34(Time, ROW) ;
HC34 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC34','f4',('Time','ROW'))
HC34[:] = TotlPoint_w_extra_12to24Z_sundy_HC34
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC34'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC34'], varattr);
        setattr(HC34, varattr, varattrVal)
        
#float HC35(Time, ROW) ;
HC35 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC35','f4',('Time','ROW'))
HC35[:] = TotlPoint_w_extra_12to24Z_sundy_HC35
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC35'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC35'], varattr);
        setattr(HC35, varattr, varattrVal)
        
#float HC36(Time, ROW) ;
HC36 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC36','f4',('Time','ROW'))
HC36[:] = TotlPoint_w_extra_12to24Z_sundy_HC36
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC36'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC36'], varattr);
        setattr(HC36, varattr, varattrVal)
        
#float HC37(Time, ROW) ;
HC37 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC37','f4',('Time','ROW'))
HC37[:] = TotlPoint_w_extra_12to24Z_sundy_HC37
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC37'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC37'], varattr);
        setattr(HC37, varattr, varattrVal)
        
#float HC38(Time, ROW) ;
HC38 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC38','f4',('Time','ROW'))
HC38[:] = TotlPoint_w_extra_12to24Z_sundy_HC38
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC38'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC38'], varattr);
        setattr(HC38, varattr, varattrVal)
        
#float HC39(Time, ROW) ;
HC39 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC39','f4',('Time','ROW'))
HC39[:] = TotlPoint_w_extra_12to24Z_sundy_HC39
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC39'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC39'], varattr);
        setattr(HC39, varattr, varattrVal)
        
#float HC40(Time, ROW) ;
HC40 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC40','f4',('Time','ROW'))
HC40[:] = TotlPoint_w_extra_12to24Z_sundy_HC40
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC40'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC40'], varattr);
        setattr(HC40, varattr, varattrVal)

#float HC41(Time, ROW) ;
HC41 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC41','f4',('Time','ROW'))
HC41[:] = TotlPoint_w_extra_12to24Z_sundy_HC41
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC41'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC41'], varattr);
        setattr(HC41, varattr, varattrVal)
        
#float HC42(Time, ROW) ;
HC42 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC42','f4',('Time','ROW'))
HC42[:] = TotlPoint_w_extra_12to24Z_sundy_HC42
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC42'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC42'], varattr);
        setattr(HC42, varattr, varattrVal)
        
#float HC43(Time, ROW) ;
HC43 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC43','f4',('Time','ROW'))
HC43[:] = TotlPoint_w_extra_12to24Z_sundy_HC43
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC43'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC43'], varattr);
        setattr(HC43, varattr, varattrVal)
        
#float HC44(Time, ROW) ;
HC44 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC44','f4',('Time','ROW'))
HC44[:] = TotlPoint_w_extra_12to24Z_sundy_HC44
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC44'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC44'], varattr);
        setattr(HC44, varattr, varattrVal)
        
#float HC45(Time, ROW) ;
HC45 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC45','f4',('Time','ROW'))
HC45[:] = TotlPoint_w_extra_12to24Z_sundy_HC45
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC45'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC45'], varattr);
        setattr(HC45, varattr, varattrVal)
        
#float HC46(Time, ROW) ;
HC46 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC46','f4',('Time','ROW'))
HC46[:] = TotlPoint_w_extra_12to24Z_sundy_HC46
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC46'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC46'], varattr);
        setattr(HC46, varattr, varattrVal)
        
#float HC47(Time, ROW) ;
HC47 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC47','f4',('Time','ROW'))
HC47[:] = TotlPoint_w_extra_12to24Z_sundy_HC47
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC47'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC47'], varattr);
        setattr(HC47, varattr, varattrVal)
        
#float HC48(Time, ROW) ;
HC48 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC48','f4',('Time','ROW'))
HC48[:] = TotlPoint_w_extra_12to24Z_sundy_HC48
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC48'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC48'], varattr);
        setattr(HC48, varattr, varattrVal)
        
#float HC49(Time, ROW) ;
HC49 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC49','f4',('Time','ROW'))
HC49[:] = TotlPoint_w_extra_12to24Z_sundy_HC49
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC49'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC49'], varattr);
        setattr(HC49, varattr, varattrVal)
        
#float HC50(Time, ROW) ;
HC50 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC50','f4',('Time','ROW'))
HC50[:] = TotlPoint_w_extra_12to24Z_sundy_HC50
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC50'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC50'], varattr);
        setattr(HC50, varattr, varattrVal)

#float HC51(Time, ROW) ;
HC51 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC51','f4',('Time','ROW'))
HC51[:] = TotlPoint_w_extra_12to24Z_sundy_HC51
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC51'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC51'], varattr);
        setattr(HC51, varattr, varattrVal)
        
#float HC52(Time, ROW) ;
HC52 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC52','f4',('Time','ROW'))
HC52[:] = TotlPoint_w_extra_12to24Z_sundy_HC52
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC52'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC52'], varattr);
        setattr(HC52, varattr, varattrVal)
        
#float HC53(Time, ROW) ;
HC53 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC53','f4',('Time','ROW'))
HC53[:] = TotlPoint_w_extra_12to24Z_sundy_HC53
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC53'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC53'], varattr);
        setattr(HC53, varattr, varattrVal)
        
#float HC54(Time, ROW) ;
HC54 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC54','f4',('Time','ROW'))
HC54[:] = TotlPoint_w_extra_12to24Z_sundy_HC54
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC54'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC54'], varattr);
        setattr(HC54, varattr, varattrVal)
        
#float HC55(Time, ROW) ;
HC55 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC55','f4',('Time','ROW'))
HC55[:] = TotlPoint_w_extra_12to24Z_sundy_HC55
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC55'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC55'], varattr);
        setattr(HC55, varattr, varattrVal)
        
#float HC56(Time, ROW) ;
HC56 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC56','f4',('Time','ROW'))
HC56[:] = TotlPoint_w_extra_12to24Z_sundy_HC56
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC56'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC56'], varattr);
        setattr(HC56, varattr, varattrVal)
        
#float HC57(Time, ROW) ;
HC57 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC57','f4',('Time','ROW'))
HC57[:] = TotlPoint_w_extra_12to24Z_sundy_HC57
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC57'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC57'], varattr);
        setattr(HC57, varattr, varattrVal)
        
#float HC58(Time, ROW) ;
HC58 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC58','f4',('Time','ROW'))
HC58[:] = TotlPoint_w_extra_12to24Z_sundy_HC58
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC58'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC58'], varattr);
        setattr(HC58, varattr, varattrVal)
        
#float HC59(Time, ROW) ;
HC59 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC59','f4',('Time','ROW'))
HC59[:] = TotlPoint_w_extra_12to24Z_sundy_HC59
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC59'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC59'], varattr);
        setattr(HC59, varattr, varattrVal)

#float HC60(Time, ROW) ;
HC60 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC60','f4',('Time','ROW'))
HC60[:] = TotlPoint_w_extra_12to24Z_sundy_HC60
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC60'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC60'], varattr);
        setattr(HC60, varattr, varattrVal)

#float HC61(Time, ROW) ;
HC61 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC61','f4',('Time','ROW'))
HC61[:] = TotlPoint_w_extra_12to24Z_sundy_HC61
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC61'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC61'], varattr);
        setattr(HC61, varattr, varattrVal)
        
#float HC62(Time, ROW) ;
HC62 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC62','f4',('Time','ROW'))
HC62[:] = TotlPoint_w_extra_12to24Z_sundy_HC62
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC62'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC62'], varattr);
        setattr(HC62, varattr, varattrVal)
        
#float HC63(Time, ROW) ;
HC63 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC63','f4',('Time','ROW'))
HC63[:] = TotlPoint_w_extra_12to24Z_sundy_HC63
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC63'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC63'], varattr);
        setattr(HC63, varattr, varattrVal)
        
#float HC64(Time, ROW) ;
HC64 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC64','f4',('Time','ROW'))
HC64[:] = TotlPoint_w_extra_12to24Z_sundy_HC64
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC64'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC64'], varattr);
        setattr(HC64, varattr, varattrVal)
        
#float HC65(Time, ROW) ;
HC65 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC65','f4',('Time','ROW'))
HC65[:] = TotlPoint_w_extra_12to24Z_sundy_HC65
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC65'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC65'], varattr);
        setattr(HC65, varattr, varattrVal)
        
#float HC66(Time, ROW) ;
HC66 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC66','f4',('Time','ROW'))
HC66[:] = TotlPoint_w_extra_12to24Z_sundy_HC66
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC66'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC66'], varattr);
        setattr(HC66, varattr, varattrVal)
        
#float HC67(Time, ROW) ;
HC67 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC67','f4',('Time','ROW'))
HC67[:] = TotlPoint_w_extra_12to24Z_sundy_HC67
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC67'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC67'], varattr);
        setattr(HC67, varattr, varattrVal)
        
#float HC68(Time, ROW) ;
HC68 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('HC68','f4',('Time','ROW'))
HC68[:] = TotlPoint_w_extra_12to24Z_sundy_HC68
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['HC68'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['HC68'], varattr);
        setattr(HC68, varattr, varattrVal)

#float PM01(Time, ROW) ;
PM01 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM01','f4',('Time','ROW'))
PM01[:] = TotlPoint_w_extra_12to24Z_sundy_PM01
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM01'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM01'], varattr);
        setattr(PM01, varattr, varattrVal)

#float PM02(Time, ROW) ;
PM02 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM02','f4',('Time','ROW'))
PM02[:] = TotlPoint_w_extra_12to24Z_sundy_PM02
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM02'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM02'], varattr);
        setattr(PM02, varattr, varattrVal)
        
#float PM03(Time, ROW) ;
PM03 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM03','f4',('Time','ROW'))
PM03[:] = TotlPoint_w_extra_12to24Z_sundy_PM03
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM03'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM03'], varattr);
        setattr(PM03, varattr, varattrVal)

#float PM04(Time, ROW) ;
PM04 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM04','f4',('Time','ROW'))
PM04[:] = TotlPoint_w_extra_12to24Z_sundy_PM04
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM04'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM04'], varattr);
        setattr(PM04, varattr, varattrVal)
        
#float PM05(Time, ROW) ;
PM05 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM05','f4',('Time','ROW'))
PM05[:] = TotlPoint_w_extra_12to24Z_sundy_PM05
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM05'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM05'], varattr);
        setattr(PM05, varattr, varattrVal)

#float PM06(Time, ROW) ;
PM06 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM06','f4',('Time','ROW'))
PM06[:] = TotlPoint_w_extra_12to24Z_sundy_PM06
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM06'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM06'], varattr);
        setattr(PM06, varattr, varattrVal)
        
#float PM07(Time, ROW) ;
PM07 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM07','f4',('Time','ROW'))
PM07[:] = TotlPoint_w_extra_12to24Z_sundy_PM07
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM07'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM07'], varattr);
        setattr(PM07, varattr, varattrVal)
        
#float PM08(Time, ROW) ;
PM08 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM08','f4',('Time','ROW'))
PM08[:] = TotlPoint_w_extra_12to24Z_sundy_PM08
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM08'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM08'], varattr);
        setattr(PM08, varattr, varattrVal)
        
#float PM09(Time, ROW) ;
PM09 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM09','f4',('Time','ROW'))
PM09[:] = TotlPoint_w_extra_12to24Z_sundy_PM09
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM09'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM09'], varattr);
        setattr(PM09, varattr, varattrVal)
        
#float PM10(Time, ROW) ;
PM10 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM10','f4',('Time','ROW'))
PM10[:] = TotlPoint_w_extra_12to24Z_sundy_PM10
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM10'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM10'], varattr);
        setattr(PM10, varattr, varattrVal)
        
#float PM11(Time, ROW) ;
PM11 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM11','f4',('Time','ROW'))
PM11[:] = TotlPoint_w_extra_12to24Z_sundy_PM11
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM11'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM11'], varattr);
        setattr(PM11, varattr, varattrVal)

#float PM12(Time, ROW) ;
PM12 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM12','f4',('Time','ROW'))
PM12[:] = TotlPoint_w_extra_12to24Z_sundy_PM12
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM12'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM12'], varattr);
        setattr(PM12, varattr, varattrVal)
        
#float PM13(Time, ROW) ;
PM13 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM13','f4',('Time','ROW'))
PM13[:] = TotlPoint_w_extra_12to24Z_sundy_PM13
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM13'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM13'], varattr);
        setattr(PM13, varattr, varattrVal)

#float PM14(Time, ROW) ;
PM14 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM14','f4',('Time','ROW'))
PM14[:] = TotlPoint_w_extra_12to24Z_sundy_PM14
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM14'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM14'], varattr);
        setattr(PM14, varattr, varattrVal)
        
#float PM15(Time, ROW) ;
PM15 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM15','f4',('Time','ROW'))
PM15[:] = TotlPoint_w_extra_12to24Z_sundy_PM15
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM15'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM15'], varattr);
        setattr(PM15, varattr, varattrVal)

#float PM16(Time, ROW) ;
PM16 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM16','f4',('Time','ROW'))
PM16[:] = TotlPoint_w_extra_12to24Z_sundy_PM16
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM16'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM16'], varattr);
        setattr(PM16, varattr, varattrVal)
        
#float PM17(Time, ROW) ;
PM17 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM17','f4',('Time','ROW'))
PM17[:] = TotlPoint_w_extra_12to24Z_sundy_PM17
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM17'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM17'], varattr);
        setattr(PM17, varattr, varattrVal)
        
#float PM18(Time, ROW) ;
PM18 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM18','f4',('Time','ROW'))
PM18[:] = TotlPoint_w_extra_12to24Z_sundy_PM18
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM18'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM18'], varattr);
        setattr(PM18, varattr, varattrVal)
        
#float PM19(Time, ROW) ;
PM19 = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('PM19','f4',('Time','ROW'))
PM19[:] = TotlPoint_w_extra_12to24Z_sundy_PM19
varattrs=["units","FieldType","MemoryOrder","description","coordinates"]
for varattr in varattrs:
    if hasattr(TotlPoint_12to24Z_sundy_file.variables['PM19'], varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file.variables['PM19'], varattr);
        setattr(PM19, varattr, varattrVal)

#char Times(Time) ;
Times = TotlPoint_w_extra_12to24Z_sundy_file.createVariable('Times','S1',('Time'))
Times[:] = TotlPoint_12to24Z_sundy_Times

#copy global attributes from TotlPoint_12to24Z_sundy_file
for varattr in TotlPoint_12to24Z_sundy_file.ncattrs():
    if hasattr(TotlPoint_12to24Z_sundy_file, varattr):
        varattrVal = getattr(TotlPoint_12to24Z_sundy_file, varattr);
        setattr(TotlPoint_w_extra_12to24Z_sundy_file, varattr, varattrVal)

TotlPoint_w_extra_12to24Z_sundy_file.close()

