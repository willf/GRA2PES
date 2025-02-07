import xarray as xr

# Directory Inputs
in_dir             = '/wrk/users/charkins/emissions/CAMSv4.2/regrid/2019_fixedOM/Month01/'
out_dir             = '/wrk/users/charkins/emissions/CAMSv4.2/regrid/2019_fixedOM/Month01/'

template_file_00z      = '/wrk/d2/charkins/CAMSv4.2/template/template_00to12Z.nc' 
template_file_12z      = '/wrk/d2/charkins/CAMSv4.2/template/template_00to12Z.nc' 


DOW                = ['weekdy','satdy','sundy']

# File inputs
cams_sector_suffix = ['_masked']
times              = ['_00z','_12z']
filetype           = '.nc4'

# Sector and Species Inputs
#species_varnames   = ['NOX','CO']
species_varnames   = ['NOX','CO','CO2','SO2','NH3','PM04','PM05','PM06','HC01',
        'HC02','HC05','HC06','HC07','HC08','HC09','HC10','HC11','HC14','HC18',
        'HC20','HC31','HC37','HC38','HC39','HC40','HC41','HC42','HC44','HC45',
        'HC48','HC49','HC50','VOC','PM25-PRI','PM10-PRI']


#species to include but have as zeros
zero_species = [ 'HC03','HC04','HC12','HC13','HC15','HC16',
         'HC17', 'HC19','HC21','HC22','HC23','HC24','HC25','HC26',
         'HC27', 'HC28','HC29','HC30','HC32','HC33','HC34',
         'HC35', 'HC36','HC43','HC46','HC47','HC51','HC52','HC53',
         'PM01', 'PM02','PM03','PM07','PM08','PM09','PM10',
         'PM11', 'PM12','PM13','PM14','PM15','PM16','PM17','PM18',
         'PM19']



sectors          = ['agriculture_livestock', 'agriculture_soils', 
                      'agriculture_wasteburning', 'power_gen', 
                      'fugitives', 'industry', 
                      'residential_commerc_comb', 'ships', 
                      'solvents', 'solidwaste_wastewater']

sectors_with_DOW = ['offroad', 'onroad'] # These are sectors that have separate DOW files


sector_fn_suffix          = {'agriculture_livestock':'_masked', 'agriculture_soils':'_masked', 
                             'agriculture_wasteburning':'_masked', 'power_gen':'_masked', 
                             'fugitives':'_masked', 'industry':'_masked', 
                             'residential_commerc_comb':'_masked', 'ships':'_unmasked', 
                             'solvents':'_masked', 'solidwaste_wastewater':'_masked',
                             'offroad':'_masked', 'onroad':'_masked'} # Lookup for suffix with masked or unmasked info

# Options
compress_output = True

# loop over the two files that make up a day 
for day_half in times:
    first = True
    # Loop over sectors without different day of week diurnals and add variables up
    for sector in sectors:
        if first:
            orig = xr.open_dataset(in_dir + sector +sector_fn_suffix[sector] +day_half + filetype, chunks={'Time': 1})
            first = False
            in_ds = orig[species_varnames]
        else:    
            in_ds = (in_ds + xr.open_dataset(in_dir + sector + sector_fn_suffix[sector] +day_half + filetype ,
                chunks={'Time': 1})[species_varnames] )
        print('Added sector: ' + sector)
        # End if
    # End first sector loop
    print('Starting DOW sector loop.')
    # Loop over each day of week
    for day in DOW:
        out_ds = in_ds.copy(deep=True)
        # Loop over sectors with different DOW profiles
        for sector in sectors_with_DOW:
            out_ds = ( out_ds + xr.open_dataset(in_dir + sector + sector_fn_suffix[sector] + '_'+ day + day_half + filetype ,
                chunks={'Time': 1})[species_varnames] )
            print('Added sector: ' + sector)
        
        # # calculate total VOC and total PM species
        # template = xr.open_dataset(template_file_00z,chunks={'Time':1})
        # out_ds['VOC']=xr.zeros_like(template['VOC'])
        # for v in mw_for_voc_total.keys():
            # out_ds['VOC'] =+ out_ds[v]*mw_for_voc_total[v]/(10**6) # convert from moles/hour to mt/hour
            
        # out_ds['PM25-PRI']=xr.zeros_like(template['PM25-PRI'])
        # out_ds['PM10-PRI']=xr.zeros_like(template['PM10-PRI'])
        # out_ds['PM25-PRI'].values = out_ds['PM04'].values + out_ds['PM05'].values + out_ds['PM06'].values # add together for totals
        # out_ds['PM10-PRI'].values = out_ds['PM04'].values + out_ds['PM05'].values + out_ds['PM06'].values # add together for totals
        
        out_ds['Times'] = xr.DataArray(data=orig['Times'].values,dims=['Time'])

        # Add in the species that have zero emissions
        if '_00z' in day_half:
            template = xr.open_dataset(template_file_00z,chunks={'Time':1})
            for species in zero_species:
                out_ds[species] = template[species]
            template.close()
        elif '_12z' in day_half: 
            template = xr.open_dataset(template_file_12z,chunks={'Time':1})
            for species in zero_species:
                out_ds[species] = template[species]
            template.close()


        # Add all of the encoding info for output
        chunk_time, chunk_lev, chunk_y, chunk_x =  1, 1, in_ds.sizes['south_north'], in_ds.sizes['west_east']
        out_coords = 'XLONG XLAT'

        if compress_output:
            data_encodingopts = {'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                                'zlib': True, 'complevel': 1, '_FillValue': None, 'coordinates': out_coords}
        else:
            data_encodingopts = {'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
                                 '_FillValue': None, 'coordinates': out_coords}
         # Set attributes based on original
        for vars in out_ds.data_vars:
            out_ds[vars].encoding = data_encodingopts
            out_ds[vars].attrs = template[vars].attrs
        out_ds.attrs = orig.attrs

                
        out_ds['XLAT'].encoding={'dtype': 'float32', '_FillValue': None}#, 'chunksizes': (834,953),'zlib': True, 'complevel': 1, '_FillValue': None}
        out_ds['XLONG'].encoding={'dtype': 'float32', '_FillValue': None}
        out_ds['Times'].encoding={'char_dim_name':'DateStrLen'}
        
        if '_12z' in day_half:
            out_fn = out_dir + day + '/' +'CAMSv4.2_12to24Z.nc'
        elif '_00z' in day_half:
            out_fn = out_dir + day + '/' +'CAMSv4.2_00to12Z.nc'
        # Save out to netCDF
        out_ds.to_netcdf(out_fn,format='netCDF4',engine='netcdf4')
        print('Saved file: ' + out_fn)
        #End second sector loop
    # End DOW loop
# End loop over each half of day

print('Program finished.')
            
            
    #         for ii in ds_out00z.data_vars:
    # ds_out00z[ii].encoding={'dtype': 'float32', 'chunksizes':(chunk_time, chunk_y, chunk_x),
    #                      'zlib': True, 'complevel': 1, '_FillValue': None}# , 'coordinates': out_coords}
            

        

            
        




