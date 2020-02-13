import numpy as np
import pdb,glob,os,sys,time,pickle
from astropy.io import fits
import matplotlib.pyplot as plt
import extract_phot 
from readcol import readcol
import querytesscut

##input file:
infile = 'targetlists/onlythegoodstuff.csv' ##put your target list here
indata = readcol(infile,fsep=',',asRecArray=True)

datadir = 'datastore/'
runsector = None ##run with runsector=None to make a first population by doing everything
runsector = None
xs,ys = 31,31 ##the tpf size to ask for, large but not crazy

##got a particular in mind? do this
##qwe  = np.where(indata.tic == 270377865)[0]
#start=qwe[0]

start = 0 ##do from the start to end of your list

##rerun things already existing?
redotpf = True
redoextract = True

##loop over targets and make lightcurves
for i in range(start,len(indata)):
    print('Up to ' + str(i+1) + ' out of ' + str(len(indata)))
    ##do a really small 1x1 cutout search to determine if something is on the chip, then do the full cutout later
    qradec_code,dddd = querytesscut.qradec(indata.tic[i],indata.ra[i],indata.dec[i],datadir=datadir,xsize=1,ysize=1,sector=None,dummymode=True)
    
    if qradec_code == -1: print('No data for this guy in this sector, skipping')
    if qradec_code ==  0: 
        if redotpf == True:
            qradec_code,tpflocations = querytesscut.qradec(indata.tic[i],indata.ra[i],indata.dec[i],datadir=datadir,xsize=xs,ysize=ys,sector=None,dummymode=False)
            if qradec_code != 0: 
                print('something very wrong here! I need your help!!!!!!!')
                pdb.set_trace()
            
        print('target observed in this sector, extracting photometry')
        tpffiles = glob.glob(tpflocations+'*_astrocut.fits')
        for j in range(0, len(tpffiles)):
            datestamp = time.strftime('%Y%m%d-%H:%M:%S', time.localtime(time.time()))
            thissector = str(int(tpffiles[j].split('-')[1].split('s')[-1]))
            ##check if we really want to re-extract this sector or not
            if runsector != None:
                sectorcheck = False
                if runsector == int(thissector): sectorcheck = True
            else: sectorcheck = True ##for the all sector run don't bother checking
            
            if sectorcheck == True:
                resultfilename  = tpflocations+'TIC'+str(indata[i].tic)+'_S'+thissector+'_FFILC'+datestamp+'.pkl'
                resultfiletrunk = tpflocations+'TIC'+str(indata[i].tic)+'_S'+thissector+'_FFILC'
                ##check for a plot directory:
                if os.path.exists(tpflocations+'plots/') == False: os.makedirs(tpflocations+'plots/')
                #remove old plots
                #pdb.set_trace()
                if  (len(glob.glob(resultfiletrunk+'*.pkl'))==0) | (redoextract == True):
                
                    
                    dfiles = glob.glob(tpflocations+'plots/TIC'+str(indata.tic[i])+'_S'+thissector+'*.pdf')
                    for df in dfiles: os.remove(df)           
                    #remove old lightcurves
                
                    lcfiles = glob.glob(tpflocations+'TIC'+str(indata.tic[i])+'_S'+thissector+'*.pkl')
                    for lc in lcfiles: os.remove(lc)
                    contig=False
                    if indata.tmag[i] <= 5.0: contig = True
                
                    ##at this point, double check the star is on the chip in this sector, and not on the very edge
                    thishdu = fits.open(tpffiles[j])
                    tot_image = thishdu[1].data[0]['FLUX']*0.0
                    for ddd in range(len(thishdu[1].data)): tot_image +=  thishdu[1].data[ddd][4]*1.0
                    zxc = np.where(np.isnan(tot_image))[0]
                
                    #pdb.set_trace()
                    if len(zxc)*1.0/len(tot_image.flatten()) < 0.05: ##if more than 5% of pixels are nans. dont even bother
                        lcdata = extract_phot.run_extraction(tpffiles[j],resultfilename,plotmovie=False,contig=contig,tmag=indata.tmag[i],sector = thissector,fixaperture=-1)
                    thishdu.close()
                    
                    ##at this point you have a lightcurve, can convert it to what I use with:
                    s3mode = False
                    if thissector == '3' : s3mode = True
                    cmdump=False
                    if indata.tmag[i] >8.0:cmdump=True
                    sectordata = extract_phot.LCconvert(lcdata,sector3=s3mode,correctmdump=cmdump)
                    outout = open(resultfilename.split('.')[0]+'_converted.pkl','wb')
                    pickle.dump(sectordata,outout)
                    outout.close()
                                        
                    #pdb.set_trace()
                else: print('Extraction Done already!, skip')
                
print('Finished Mate')
