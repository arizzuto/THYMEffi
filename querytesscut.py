import pdb,os,time,subprocess,zipfile

##a function to download tesscut TPF's

def qtic(ticnum,datadir = '/Volumes/UTOld/tessdata/',xsize=31,ysize=31,sector=None,dummymode=False):
    import tic_query
    ##get the ra and dec from TIC then just run the other function
    ra,dec = tic_query.tic_radec(ticnum)
    dummy,dummy2 = qradec(ticnum,ra,dec,datadir=datadir,xsize=xsize,ysize=ysize,sector=sector,dummymode=dummymode)
    return dummy,dummy2
    
    
def qradec(tic,ra,dec,datadir = '/Volumes/UTOld/tessdata/',xsize=31,ysize=31,sector=None,dummymode=False):
    ##ra and dec are either floats or strings
    ##sector and tic are integeres, as as x/ysize
    
    ##build the command
    queryurl = 'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra='+str(ra)+'&dec='+str(dec)+'&y='+str(ysize)+'&x='+str(xsize)
    if sector != None: queryurl += '&sector='+str(sector)
    command  = '/usr/local/bin/wget -O ./tesscut_tmp/latest.zip "'+queryurl + '"'  
    
    process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)    
    process.wait()
    #print('subprocess returncode: ' + str(process.returncode))
    
    if process.returncode !=0:
        ##try the command again
        process = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)    
        process.wait()
    
        if process.returncode != 0: 
            print('bad return code, something wrong!')
            errorfile = open('latest_tesscut_errorfile.txt','wb')
            errorfile.write(command + ' \n')
            errorfile.write(str(process.returncode) + ' \n')
            errorfile.flush()
            errorfile.close()
            
            pdb.set_trace()
        
        
        
    #pdb.set_trace()
    ##'/usr/local/bin/wget -O ./tesscut_tmp/latest.zip "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=179.2008&dec=-22.4893&y=1&x=1&sector=10"'
    ##'/usr/local/bin/wget -O ./tesscut_tmp/latest.zip "https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=347.45642&dec=-14.51055&y=1&x=1&sector=10"'
    filesize = os.path.getsize('./tesscut_tmp/latest.zip')
    if filesize < 5000:
        ##this is the case for a non-observed target in the input sector
        return -1,-1
    if dummymode == True: return 0,0 ##if in dummymode=True, just checking target observation, dont extract any files
    
    ##now repackage the download, basically we want one subdir for each tic number that gets updated every time if new stuff appears   
    ##check directory exists:
    datalocation = datadir+'tic'+str(tic)+'/'
    if os.path.exists(datalocation)==False: os.makedirs(datalocation)
    
    ##unzip to the new data location
    zip_ref = zipfile.ZipFile('./tesscut_tmp/latest.zip', 'r')
    zip_ref.extractall(datalocation)
    zip_ref.close()
    ##now write in the corresponding log
    datestamp = time.strftime('%Y%m%d-%H:%M:%S', time.localtime(time.time()))
    logfile = open(datalocation+'tic'+str(tic)+'_dllog.txt','ab')
    logfile.write(datestamp+' '+str(ra)+' '+str(dec)+' '+str(xsize)+'x'+str(ysize)+' '+str(sector)+' \n')
    logfile.close()

    return 0,datalocation
    
    
    
##uncomment to test with DS tuc
# tic = 410214986
# ra = 354.91395712
# dec = -69.19558694
# test = qradec(tic,ra,dec)

##uncomment to test an unobserved target:
# tic = 17554529
# ra = 67.15464437
# dec=19.18028285
# test = qradec(tic,ra,dec)