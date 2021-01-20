#!/usr/bin/env python
#******************************************************************************
#  Name:     s2cloudless.py (colab version)
#  Purpose:  Create a cloud-free image from GEE
#  Usage (from command line):             
#    python s2cloudless.py  [options] areaOfInterest outFilename
#
#  Copyright (c) 2021 Mort Canty

import ee

ee.Initialize()

    
def main(): 

    usage = '''            
Usage: 
--------------------------------------

Deep learning classification of S2 images

python %s [OPTIONS] aoi outfilename

Options:
  -h            this help
  -s  <string>  start date (default '2020-06-01')
  -e  <string>  end date (default '2020-09-01')
  
  -------------------------------------'''%sys.argv[0]            
                    
    options,args = getopt.getopt(sys.argv[1:],'hnm:t:d:s:v:')

    start_date = '2020-06-01'
    end_date = '2020-09-01'
    CLOUD_FILTER = 60
    CLD_PRB_THRESH = 40
    NIR_DRK_THRESH = 0.15
    CLD_PRJ_DIST = 2
    BUFFER = 100
    
    for option, value in options: 
        if option == '-h':
            print(usage)
            return 
        elif option == '-s':
            start_date = value
        elif option == '-e':
            end_date = value) 
    if len(args) != 2:
        print( 'Incorrect number of arguments' )
        print(usage)
        sys.exit(1)  
        
    outfile = sys.
        
    s2_sr_cld_col = get_s2_sr_cld_col(aoi, start_date, end_date)        
        

if __name__ == '__main__':
    main()    
    
    