import numpy as np
import sys
import os
from extract_counts import *
from BEHR_countbins import *
from run_BEHR import *
from re import sub
from create_behr_files import region_area

def HR_main(obsid,position,N,lines,divide_energy=2000,override=False,subtract_start=True):
    #BEHR_DIR = '/data/reu/slarner/BEHR_contain/BEHR'
    #or use:
    BEHR_DIR = "/home/ievajankute/Documents/Work/BEHR"

    download_obsid(obsid)

    primary = f'./{obsid}/primary'

    #remove the colons from the position
    position_basic = sub('\:','',position)
    working_dir = f'{primary}/{position_basic}'

    if override:
        src_region = input('Enter path to src region file: ').strip()
        bkg_region = input('Enter path to bkg region file: ').strip()

    else:
        print('Running srcflux to make regions...')
        #If not overriden, we make the regions with src flux
        make_regions(obsid,position,f'{working_dir}/')

        src_region = unglob(glob.glob(f'{working_dir}/*srcreg.fits'),True)
        print(src_region)
        bkg_region = unglob(glob.glob(f'{working_dir}/*bkgreg.fits'),True)

    evt = unglob(glob.glob(f'{primary}/*evt2*'))
    
    #check to see if the background is wrong
    #when the srcflux source region extends over the edge of a chip, the
    #corresponding background region will be far far too large
    #I am going to detect this by comparing the area of each

    try:

        area = region_area(evt,bkg_region,1000)/region_area(evt,src_region,1000)
    except Exception as e:
        print('If OSError because the virtual file for a region could not be opened')
        print('Remember to save the bkg region in ciao format instead of ds9')
        print('Talk me if you need help')

        raise e

    print('Making BEHR bash file...')

    outfile = f'{working_dir}/BEHR_bash.txt'

    BEHR_outdir = f'{BEHR_DIR}/{obsid}/{position_basic}'

    subprocess.run(f'rm -rf {BEHR_outdir}',shell=True)
    os.makedirs(BEHR_outdir)

    make_behr(evt,src_region,bkg_region,divide_energy,BEHR_DIR,outfile,BEHR_outdir,N)

    print('Running BEHR...')

    run_BEHR(outfile)

    print('Making plot...')

    if subtract_start:
        dmextract.punlearn()
        dmextract.infile = f'{evt}[bin time=::3.24014]'
        dmextract.clobber = 'yes'
        dmextract.opt = 'generic'
        dmextract.outfile = 'temp.fits'
        dmextract()

        dmlist.punlearn()
        dmlist.infile = 'temp.fits[cols time]'
        dmlist.opt = 'data, clean'

        start_time = float(dmlist().splitlines()[1])
    else:
        start_time = 0

    time,uppers,lowers,meds = plot_BEHR_constcounts(BEHR_outdir,N,position_basic,obsid,evt,src_region,lines=lines,start_time=start_time,show=True)
    
    """
    print('Saving...')
    to_save = np.column_stack((time,meds,uppers,lowers))

    np.savetxt(f'./HR_saved_{position}_{obsid}.txt',to_save,delimiter=',',header='Time,Median,Upper Error,Lower Error')

    return to_save #Change to return plotting parameters
    """
    return time,uppers,lowers,meds