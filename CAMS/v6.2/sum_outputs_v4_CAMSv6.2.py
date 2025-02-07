import xarray as xr

class inputs():

    # Directory Inputs
    in_dir             = '/wrk/users/charkins/emissions/CAMSv6.2/regrid_nomask/[YYYY]/Month[MM]/'
    out_dir             = '/wrk/users/charkins/emissions/CAMSv6.2/regrid_nomask/[YYYY]/Month[MM]/'

    template_file_00z      = '/wrk/d2/charkins/templates/with_additional_VCPs_ffCO2/template_new_00to12Z.nc'
    template_file_12z      = '/wrk/d2/charkins/templates/with_additional_VCPs_ffCO2/template_new_12to24Z.nc'


    DOW                = ['weekdy','satdy','sundy']

    # File inputs
    times              = ['_00z','_12z']
    filetype           = '.nc4'

    # Sector and Species Inputs
    #species_varnames   = ['NOX','CO']
    species_varnames   = ['NOX','CO','ffCO2','CO2','SO2','NH3','PM04','PM05','PM06','HC01',
            'HC02','HC05','HC06','HC07','HC08','HC09','HC10','HC11','HC14','HC18',
            'HC20','HC31','HC37','HC38','HC39','HC40','HC41','HC42','HC44','HC45',
            'HC48','HC49','HC50']

    #species to include but have as zeros
    zero_species = [ 'HC03','HC04','HC12','HC13','HC15','HC16',
             'HC17', 'HC19','HC21','HC22','HC23','HC24','HC25','HC26',
             'HC27', 'HC28','HC29','HC30','HC32','HC33','HC34',
             'HC35', 'HC36','HC43','HC46','HC47','HC51','HC52','HC53',
             'HC54','HC55','HC56','HC57','HC58','HC59','HC60','HC61',
             'HC62','HC63','HC64','HC65','HC66','HC67','HC68',
             'PM01', 'PM02','PM03','PM07','PM08','PM09','PM10',
             'PM11', 'PM12','PM13','PM14','PM15','PM16','PM17','PM18',
             'PM19','PM10-PRI','PM25-PRI']



    sectors          = ['agriculture_livestock', 'agriculture_soils', 
                          'agriculture_wasteburning','commerc_comb', 'power_gen', 
                          'fugitives', 'industry', 
                          'refineries','residential_comb', 'ships', 
                          'solvents', 'solidwaste_wastewater','solidwaste_landfill']

    sectors_with_DOW = ['offroad', 'onroad'] # These are sectors that have separate DOW files


    sector_fn_suffix          = {'agriculture_livestock':'_unmasked', 'agriculture_soils':'_unmasked', 
                                 'agriculture_wasteburning':'_unmasked','commerc_comb':'_unmasked', 'power_gen':'_unmasked', 
                                 'fugitives':'_unmasked', 'industry':'_unmasked', 
                                 'residential_comb':'_unmasked','refineries':'_unmasked', 'ships':'_unmasked', 
                                 'solvents':'_unmasked', 'solidwaste_wastewater':'_unmasked','solidwaste_landfill':'_unmasked',
                                 'offroad':'_unmasked', 'onroad':'_unmasked'} # Lookup for suffix with masked or unmasked info
    
    #years= [2019]
    years = [2023] #2006,2007,2008,2009,2011,2012,2013,2014,2016,2017,2018
    months = [1,2,3,4,5,6,7,8,9,10,11,12]
    
    # Options
    compress_output = True

def make_files(year,month):
    # loop over the two files that make up a day 
    for day_half in inputs.times:
        first = True
        # Loop over sectors without different day of week diurnals and add variables up
        for sector in inputs.sectors:
            if first:
                orig = xr.open_dataset(inputs.in_dir.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)) + sector +inputs.sector_fn_suffix[sector] +day_half + inputs.filetype, chunks={'Time': 1},cache=False)
                first = False
                in_ds = orig[inputs.species_varnames]
            else:    
                in_ds = (in_ds + xr.open_dataset(inputs.in_dir.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)) + sector + inputs.sector_fn_suffix[sector] +day_half + inputs.filetype ,
                    chunks={'Time': 1},cache=False)[inputs.species_varnames] )
            print('Added sector: ' + sector)
            # End if
        # End first sector loop
        print('Starting DOW sector loop.')
        # Loop over each day of week
        for day in inputs.DOW:
            out_ds = in_ds.copy(deep=True)
            # Loop over sectors with different DOW profiles
            for sector in inputs.sectors_with_DOW:
                out_ds = ( out_ds + xr.open_dataset(inputs.in_dir.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)) + sector + inputs.sector_fn_suffix[sector] + '_'+ day + day_half + inputs.filetype ,
                    chunks={'Time': 1})[inputs.species_varnames] )
                print('Added sector: ' + sector)
            
            out_ds['Times'] = xr.DataArray(data=orig['Times'].values,dims=['Time'])

            # Add in the species that have zero emissions
            if '_00z' in day_half:
                template = xr.open_dataset(inputs.template_file_00z,chunks={'Time':1})
                for species in inputs.zero_species:
                    out_ds[species] = template[species]
                template.close()
            elif '_12z' in day_half: 
                template = xr.open_dataset(inputs.template_file_12z,chunks={'Time':1})
                for species in inputs.zero_species:
                    out_ds[species] = template[species]
                template.close()


            # Add all of the encoding info for output
            chunk_time, chunk_lev, chunk_y, chunk_x =  1, 1, in_ds.sizes['south_north'], in_ds.sizes['west_east']
            out_coords = 'XLONG XLAT'

            if inputs.compress_output:
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
                out_fn = inputs.out_dir.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)) + day + '/' +'CAMSv6.2_12to24Z.nc'
            elif '_00z' in day_half:
                out_fn = inputs.out_dir.replace('[YYYY]',str(year)).replace('[MM]',str(month).zfill(2)) + day + '/' +'CAMSv6.2_00to12Z.nc'
            # Save out to netCDF
            out_ds.to_netcdf(out_fn,format='netCDF4',engine='netcdf4')
            print('Saved file: ' + out_fn)
            #End second sector loop
        # End DOW loop
    # End loop over each half of day
def main():

    for year in inputs.years:
        for month in inputs.months:
            print('Working on emissions for: '+ str(year) + str(month).zfill(2))
            
            make_files(year,month)
        
    print('Program finished.')
            
if __name__ == "__main__":
    main()          

            

        

            
        




