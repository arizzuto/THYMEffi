import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pdb,os,sys,time,pickle
import mpyfit
sys.path.append('../K2pipe')
##import k2sff ##not needed anymore beyond testing
from astropy.io import fits
from skimage.feature import peak_local_max
from scipy import interpolate

def robustmean(y,cut):     
    absdev    = np.absolute(y-np.median(y))
    medabsdev = np.median(absdev)/0.6745
    if medabsdev < 1.0e-24: medabsdev = np.mean(absdev)/0.8
    gind = np.where(absdev <= cut*medabsdev)[0]
    sigma    = np.std(y[gind])
    sc       = np.max([cut,1.0])
    
    if sc <= 4.5: sigma = sigma/(-0.15405+0.90723*sc-0.23584*sc**2+0.020142*sc**3) 
    goodind  = np.where(absdev <= cut*sigma)[0]
    #goodpts  = y[goodind]
    sigma    = np.std(y[goodind])
    
    sc = np.max([cut,1.0])
    if sc <= 4.5: sigma=sigma/(-0.15405+0.90723*sc-0.23584*sc**2+0.020142*sc**3)      
    return np.mean(y[goodind]),sigma/np.sqrt(len(y)-1.0),len(y) - len(goodind),goodind
    
def lcbin(time2,lc,nbins,usemean=False,userobustmean=False,linfit=False):      

    time = time2-np.nanmin(time2)
    
    usemedian = (1-usemean) & (1-userobustmean)
    
    binsize = (np.nanmax(time)-np.nanmin(time))/float(nbins)
    fbin    = np.empty(nbins)*np.nan
    ebin    = np.empty(nbins)*np.nan
    tbin    = (np.arange(nbins,dtype=float)+0.5)*binsize + np.nanmin(time2)
    if  (userobustmean == True) or (linfit==True): allgoodind = np.array([-1])


    for i in range(nbins):
        w = np.where((time >= i*binsize) & (time < (i+1)*binsize))[0]
        if len(w) > 0:
            if usemedian == True:
                fbin[i] = np.nanmedian(lc[w])
                if len(w) > 1: ebin[i] = 1.253*np.std(lc[w])/np.sqrt(len(w))
            
            if usemean == True:
                fbin[i] = np.nanmean(lc[w])
                if len(w) > 1: ebin[i] = np.std(lc[w])/np.sqrt(len(w))
            
            if userobustmean == True:
                fbin[i], j1,j2,goodind = robustmean(lc[w],3)
                if len(w) > 1:
                    if len(goodind) > 1: 
                        ebin[i] = np.std(lc[w[goodind]])/np.sqrt(len(w))
                    else: ebin[i] = np.inf
                allgoodind = np.concatenate((allgoodind,w[goodind]))
        
            if linfit == True:
                if len(w) > 2:
                    bg,rms,goodind,polyfit = robustpolyfit(time2[w],lc[w],1)
                    fbin[i] = np.polyval(polyc,tbin[i])
                    ebin[i] = rms/np.sqrt(len(w))
                    allgoodind = np.concatenate((allgoodind,w[goodind]))
                if len(w) <= 1: ##if not enough points revert to robust mean
                    fbin[i],j1,j2,goodind = robustmean(lc[w],3,goodind)
                    if len(w) >1:  ##this doesnt make sense!!!!
                        ebin[i] = np.std(lc[w[goodind]])/np.sqrt(len(w))
                        allgoodind = np.concatenate((allgoodind,w[goodind]))
                
    
    if robustmean == True: allgoodind = allgoodind[1:]
    return fbin,binsize,ebin,allgoodind,tbin
    
def clipfitline(t,f,jumpspot,niter=5,nojump=False):
    qwe = np.where(t>-np.inf)[0]
    
    line = np.polyfit(t[qwe],f[qwe],2)

    left = np.where(t<jumpspot)[0]
    right= np.where(t>jumpspot)[0]
    #pdb.set_trace()
    ssssss=np.nanmedian(f[right])-np.nanmedian(f[left])
    pstart = np.array([line[0],line[1],line[2],0.001,jumpspot])
    lerr = np.nanmedian(np.absolute(f[qwe]-np.polyval(line,t[qwe])))
    smod,sres=mdump_line(pstart,(t[qwe],f[qwe],f[qwe]*0.0+lerr),model=True)

    parinfo2 = [{'fixed':False, 'limits':(None,None),'step':0.001} for dddd in range(pstart.shape[0])]
    parinfo2[4]['fixed']=True
    if nojump == True:
        pstart[3]=0.000
        parinfo2[3]['fixed']=True
        
    for j in range(5):
        smod,sres=mdump_line(pstart,(t[qwe],f[qwe],f[qwe]*0.0+lerr),model=True)

        llpars,llext=mpyfit.fit(mdump_mod,pstart,args=(t[qwe],f[qwe],f[qwe]*0.0+lerr),parinfo=parinfo2)
        lmod,lres=mdump_mod(llpars,(t[qwe],f[qwe],f[qwe]*0.0+lerr),model=True)
        lmad = np.nanmedian(np.absolute(f[qwe]-lmod))
        keep = np.where((np.absolute(lmod-f[qwe])/lmad<10.0))[0]
        if j < niter-1:qwe=qwe[keep]
        pstart = llpars
        lerr=lmad
    
    ##give the full model:
    ##recalculate lres
    #lmod,lres=mdump_mod(llpars,(t[qwe],f[qwe],f[qwe]*0.0+lerr),model=True)
    themodel,dummy = mdump_mod(llpars,(t,f,f*0.0+lerr),model=True)
    #pdb.set_trace()
    return t,themodel,lres,qwe,llpars,lerr

def mdump_mod(p,args,model=False):
    #p=(a2,a1,a0,jump,jumpspot)
    t,flx,sigf=args
    themodel=t**2*p[0]+t*p[1]+p[2]
    xxx = np.where(t>p[4])
    themodel[xxx]+=p[3]
    if model == True: return themodel,(flx-themodel)/sigf
    return (flx-themodel)/sigf


def mdump_line(p,args,model=False):
    #p=(a2,a1,a0,jump)
    t,flx,sigf=args
    themodel=t**2*p[0]+t*p[1]+p[2]
    
    if model == True: return themodel,(flx-themodel)/sigf
    return (flx-themodel)/sigf
    




def gaussmodel_old(p,args,model=False):
    '''
    simple 2D gaussian model with background
    ##params are p=[x0,y0,sigx,sigy,amp,bg,crossterm]
    '''
    imx,imy,im,sig_im = args
    #pdb.set_trace()
    gaussmod = np.exp(-(imx-p[0])*(imy-p[1])*p[6])*np.exp(-(imx-p[0])**2/2/p[2]**2)*np.exp(-(imy-p[1])**2/2/p[3]**2)*p[4] + p[5]
    if np.isnan(gaussmod).any() == True:
        pdb.set_trace()    
    resid = (im-gaussmod)/sig_im
   # print 'calling'
    #pdb.set_trace()
    if model == True: 
        return gaussmod, resid
    else: return resid.flatten()

def gaussmodel(p,args,model=False):
    '''
    simple 2D gaussian model with background
    ##params are p=[x0,y0,sigx,sigy,amp,bg,crossterm]
    '''
    imx,imy,im,sig_im = args
    #pdb.set_trace()
    ##gaussmod = np.exp(-(imx-p[0])*(imy-p[1])*p[6]/(1-p[6]**2)/p[2]/p[3])*np.exp(-(imx-p[0])**2/2/p[2]**2/(1-p[6]**2))*np.exp(-(imy-p[1])**2/2/p[3]**2/(1-p[6]**2))*p[4] + p[5]
    gaussmod = np.exp(-(imx-p[0])*(imy-p[1])*p[6]/(1-p[6]**2)/p[2]/p[3] - (imx-p[0])**2/2/p[2]**2/(1-p[6]**2) - (imy-p[1])**2/2/p[3]**2/(1-p[6]**2))*p[4]+p[5]
    
    if np.isnan(gaussmod).any() == True:
        pdb.set_trace()    
    resid = (im-gaussmod)/sig_im
   # print 'calling'
    #pdb.set_trace()
    if model == True: 
        return gaussmod, resid
    else: return resid.flatten()



def tess_gauss_fit(im,sig_im,centroidmode='fixed',xsigcmode=0.5,stopper=False):

    ##make a cutout
    ##evertythin more than 5 pixels away should be ignored in fit
    

    ##first make the x/y grids
    xvect = np.arange(im.shape[0])
    yvect = np.arange(im.shape[1])
    
    imx,imy = np.meshgrid(yvect,xvect)
    #pdb.set_trace()
    if len(xvect) % 2 != 0:
        xstart = len(xvect)/2 
    else: xstart = len(xvect)/2-1
    
    
    
    #
    ##adjust the starting position to the brightest pixel in a small distance from the initial guess
    shiny = np.unravel_index(np.argmax(im[xstart-2:xstart+2,xstart-2:xstart+2]),im[xstart-2:xstart+2,xstart-2:xstart+2].shape)
    starterx,startery = xstart-2+shiny[0],xstart-2+shiny[1]
    #pdb.set_trace()
    #if centroidmode=='free':
    #    xxx =np.where(im == np.max(im))
    #    ystart = xxx[0][0]
    #    xstart = xxx[1][0]
        
    ##pdb.set_trace()
    qwe=  np.where((imx-starterx)**2 + (imy-startery)**2 > 4**2)
    qwe2 =  np.where((imx-starterx)**2 + (imy-startery)**2 <= 2**2)
    sig_im[qwe[0],qwe[1]] = im[qwe[0],qwe[1]]*10.0
    pstart = np.array([imx[starterx,startery]*1.0,imy[starterx,startery]*1.0,0.5,0.5,np.nanmax(im[qwe2[0],qwe2[1]]),np.max([np.nanmin(im[qwe2[0],qwe2[1]]),0.1]),0.01])

    qwe = np.where(np.isnan(im))
    im[qwe[0],qwe[1]] = 0.0
    sig_im[qwe[0],qwe[1]] = np.nanmax(im)*100
   

    ##limit certain parameters
    parinfo2 = [{'fixed':False, 'limits':(None,None)} for dddd in range(pstart.shape[0])]
    parinfo2[4]['limits']    = (0.01,None)
    parinfo2[5]['limits']    = (0.01,None)
    parinfo2[6]['limits']    = (0.0001,0.99)
    parinfo2[0]['limits']    = (pstart[0]-.2,pstart[0]+.2)
    parinfo2[1]['limits']    = (pstart[1]-.2,pstart[1]+.2)
    parinfo2[2]['limits']    = (0.1,1.0)
    parinfo2[3]['limits']    = (0.1,1.0)
    if centroidmode == 'free': ##let the centroid move alot more, mostly for custom stuff
    
        parinfo2[0]['limits']    = (pstart[0]-1.0,pstart[0]+1.0)
        parinfo2[1]['limits']    = (pstart[1]-1.0,pstart[1]+1.0)
        parinfo2[2]['fixed'] = True
        parinfo2[3]['fixed'] = True
       # parinfo2[5]['fixed'] = True
        ##parinfo[6]['fixed'] = True
        pstart[2]  = xsigcmode
        pstart[3]  = xsigcmode
        pstart[5] = 10.0
        #pstart[0] = starterx
        #pstart[1] = startery
        
    #pdb.set_trace()  
    ##parinfo2[6]['fixed'] = True
   
    if stopper == True: pdb.set_trace()
    try:
        fitpars,extrares = mpyfit.fit(gaussmodel,pstart,args=(imx,imy,im,sig_im),maxiter=200, parinfo=parinfo2)
        bestmod,bestresid = gaussmodel(fitpars,(imx,imy,im,sig_im),model=True)
    except: 
        print('the image fit failed!')
        pdb.set_trace()
    return fitpars[0],fitpars[1],fitpars[5],fitpars[2],fitpars[3]

def tess_circmask(im,sig_im,r,xcen,ycen):
    ##first make the x/y grids
    xvect = np.arange(im.shape[0])
    yvect = np.arange(im.shape[1])
    imx,imy = np.meshgrid(yvect,xvect)
   # pdb.set_trace()
    inradius = np.where(np.sqrt((imx-xcen)**2 + (imy-ycen)**2) <= r)
    thisflux = np.sum(im[inradius[0],inradius[1]])/len(inradius[0])
    thiserror = np.sqrt(np.sum(sig_im[inradius[0],inradius[1]]**2))/len(inradius[0])
    return thisflux,inradius[0],inradius[1],imx,imy,thiserror

def cstep(map,shape,xstart,ystart):

    shape[xstart,ystart] = 1
    xm = map.shape[0]-1
    ym = map.shape[1]-1
    #pdb.set_trace()
    if (xstart -1 > 0)  & (map[xstart-1,ystart] >0) & (shape[xstart-1,ystart] == 0): shape = cstep(map,shape,xstart-1,ystart)
    if (xstart +1 < xm) & (map[xstart+1,ystart] >0) & (shape[xstart+1,ystart] == 0): shape = cstep(map,shape,xstart+1,ystart)
    
    if (ystart -1 > 0)  & (map[xstart,ystart-1] >0) & (shape[xstart,ystart-1] == 0): shape = cstep(map,shape,xstart,ystart-1)
    if (ystart +1 < ym) & (map[xstart,ystart+1] >0) & (shape[xstart,ystart+1] == 0): shape = cstep(map,shape,xstart,ystart+1)
    return shape

def contiguousregions(image,xstart,ystart,level=2,ylim =8): ####note that level=1 seemed ok, 5 had some jumps
    ##code to find contiguous regions of pixels, useful for saturated images.
    ##like AV's code
    nimage = (image - np.nanmedian(image)) / (np.nanmax(image)-np.nanmedian(image))
    sigma = np.std(nimage)/len(nimage.flatten())*1.43
    hi = np.where(nimage >= np.nanmedian(nimage)+level*sigma)
    #lo = np.where(nimage < np.median(nimage)+level*sigma)
    mapimage = nimage*0.0
    mapimage[hi[0],hi[1]] = 1.0
    shape = mapimage*0.0
    shape = shape.astype(int)
    cshape = cstep(mapimage,shape,xstart,ystart)
    cshape[:,0:np.max([ystart-ylim,0])] = 0
    cshape[:,np.min([ystart+ylim,cshape.shape[1]]):]  = 0
    return cshape
    
def run_extraction(filename,resultfilename,aperture_radius=3,plotmovie=False,contig=False,tmag=10.,sector=1,centroidmode='fixed',fixaperture=-1,bgtype='percentile',noplots=False,nosave=False,return_image=False):
    '''
    give this the tpf filename location and and output filename location to save the data to
    the aperture is up to the user,     
    '''
    if (bgtype != 'linear') & (bgtype != 'cubic') & (bgtype != 'percentile') & (bgtype != 'simple'): 
        print('select a valid backgrounding method: | linear cubic percentile(default) simple | and retry')
        pdb.set_trace()
    
    hdu      = fits.open(filename)
    xcen     = np.zeros(len(hdu[1].data))
    ycen     = np.zeros(len(hdu[1].data))
    bg       = np.zeros(len(hdu[1].data))
    medim    = np.zeros((len(hdu[1].data),20))
    standard_medim = np.zeros(len(hdu[1].data))
    interpbg = np.zeros(len(hdu[1].data))
    fluxcirc = np.zeros(len(hdu[1].data))
    fluxcirc_err = np.zeros(len(hdu[1].data))

    qual = np.zeros(len(hdu[1].data),dtype=int)
    head = hdu[1].header

    cadtime  = np.zeros(len(hdu[1].data))
    tot_image = hdu[1].data[0][4]*0.0
    tot_image_err = hdu[1].data[0][5]*0.0
    for i in range(len(hdu[1].data)): 
        qual[i]    = hdu[1].data[i]['QUALITY']*1
        
        ##this qualifier: (np.isnan(hdu[1].data[i][4]).any() == False) & not needed I think
        if  (qual[i] == 0) & (np.nanmax( hdu[1].data[i][4]) > 0.0): 
            tot_image +=  hdu[1].data[i][4]*1.0
            tot_image_err +=  hdu[1].data[i][5]*1.0
        nannan =np.where(np.isnan(tot_image) == True)
        tot_image[nannan[0],nannan[1]] = np.nanmedian(tot_image)
        
    #find all stars
    
    threshold_mult=0.01

    qwe =np.where(tot_image > np.nanmax(tot_image)*threshold_mult)
    if len(qwe[0]) > 0.9*len(tot_image.flatten()): 
        badmult=True
        while badmult == True:     
            threshold_mult += 0.005
            qwe =np.where(tot_image > np.nanmax(tot_image)*threshold_mult)
            if len(qwe[0]) < 0.9*len(tot_image.flatten()): badmult=False
    scord = peak_local_max(tot_image, min_distance=1,exclude_border=False,threshold_rel=threshold_mult,num_peaks=10) ##=np.median(tot_image)+np.std(tot_image)*1.25/np.sqrt(len(tot_image.flatten()))*3
    
    ##now create a list of coords that should be ignored in backgrounding
    starcoords = ()
    windsize = 6
    for cc in scord:
        for c1 in range(windsize*2+1):
            for c2 in range(windsize*2+1):
              #  pdb.set_trace()
                thiscoord = [cc[0]-windsize+c1,cc[1]-windsize+c2]
                if (thiscoord[0] >= 0) & (thiscoord[1]>=0) & (thiscoord[0] < tot_image.shape[0]) & (thiscoord[1]<tot_image.shape[1]):  
                    if tot_image[thiscoord[0],thiscoord[1]] >= np.nanmin([np.nanmax(tot_image)*0.01,np.median(tot_image)*1.1]): ###only mask out pixel aboce the threshold level., dont over extend windows
                        starcoords += (thiscoord,)
                  

    starcoords = np.array(starcoords)

    ##for i in range(len(hdu[1].data)): tot_image_err +=  hdu[1].data[i][5]*1.0
    ##if the image was cut from the edge of a chip, correct zeros otherwise algorithms can fail
    #pdb.set_trace()
    if (np.nanmax(tot_image) == 0) | (len(np.where(np.isnan(tot_image.flatten())==True)[0]) == len(tot_image.flatten())):
        print( 'bad image, not data')
        pdb.set_trace()
        return -1
    edged = np.where(tot_image == 0)
    other = np.where(tot_image.flatten() >  0)
    #tot_image[edged[0],edged[1]] = np.nanmedian(tot_image.flatten()[other])
    tot_image_err[edged[0],edged[1]] = np.nanmedian(tot_image_err.flatten()[other])*100 ##large errors on the missing pixels

    if contig == True:
        ##now find ajdacent pixels that are above some value
       ## pdb.set_trace()
        cshape = contiguousregions(tot_image,tot_image.shape[0]/2+1,tot_image.shape[0]/2+1)
        zxc = np.where((cshape !=1) & (tot_image>0)) ##again remove background points that would be zeros
        rshape = cshape*0.0
        rshape[zxc[0],zxc[1]] = 1

    if contig == False:
        
        xcen1,ycen1,fakebg,xsig,ysig = tess_gauss_fit(tot_image*1.0,tot_image_err*1.0,centroidmode=centroidmode)
       # pdb.set_trace()
        ##failure mode: target is way too faint for tess and the gaussian size is tiny--> no pixels included:
        ##set to 2 pixels min
        #pdb.set_trace()
        if (xsig<0.5) | (ysig<0.5): ysig,xsig = 0.5,0.5  
        if xsig > 0.7 : xsig=0.7
        if tmag > 12 :
            if xsig > 0.6: xsig = 0.6  
        
        ##if the aperture size is fixed, ignore the numbers produced above:
        if fixaperture >0:
            xsig = fixaperture/5.0

    for i in range(len(hdu[1].data)):
        #print str(i) + ' of ' + str(len(hdu[1].data))
        if plotmovie == True :plt.clf()
        im = hdu[1].data[i][4]
        im_err = hdu[1].data[i][5]

        
        cadtime[i] = hdu[1].data[i][0]*1.0           
        if contig == False:
            xcen[i] = xcen1
            ycen[i] = ycen1
            fluxcirc[i],xin,yin,xim,yim,fluxcirc_err[i] = tess_circmask(im,im_err,np.round(xsig*5),xcen1,ycen1)  
            if centroidmode == 'free': ##return the centroid positions for each image if centroid mode is set to free, for custom stuff really
                if np.isnan(im).any() == False:
                    stopper=False
                    #if i == 82: stopper=True
                    thisxcen,thisycen,fakebg,xsigdummy,ysigdummy = tess_gauss_fit(im,im_err,centroidmode=centroidmode,xsigcmode=xsig,stopper=stopper)
                    xcen[i]=thisxcen
                    ycen[i]=thisycen
                    print('image ' + str(i+1) + ' out of ' + str(len(hdu[1].data)))
                fluxcirc[i],xin,yin,xim,yim,fluxcirc_err[i] = tess_circmask(im,im_err,np.round(xsig*5),xcen1,ycen1)  
                      
            
            cshape = im*0
            cshape[xin,yin] = 1
            cshape = cshape.astype(int)
            
            ##find the lowest 5% of pixels:
            ##if there are zeros in the image, need to ignore those for percentile to make sense
             ##we should also make all the pixels in the aperture zero with a boundary a pixel wider than the size of the aperture
            ddd,xin2,yin2,xim2,yim2,dddd = tess_circmask(im,im_err,np.round(xsig*5)*1+1,xcen1,ycen1)  
            imim = im*1.0
            imim[xin2,yin2] = 0.0
            #pdb.set_trace()
            flatim = imim.flatten()
            nonzero = np.where(flatim >0)[0]
            imim2 = im*1.0
            imim2[starcoords[:,0],starcoords[:,1]] = np.nan
            imim2[xin2,yin2] = np.nan
            standard_medim[i] = np.nanmean(imim2)
            okcoords = np.array(np.where(imim > 0.0)).T
            okcoords2 = np.array(np.where((np.isnan(imim2) == False) & (imim>0.0)) ).T
            xvector = np.arange(im.shape[0])
            yvector = np.arange(im.shape[1])
            imxg,imyg = np.meshgrid(xvector,yvector)
           # if i == 15: pdb.set_trace()
         #   pdb.set_trace()
            if (bgtype == 'linear') | (bgtype == 'cubic'):  
                if np.nanmax(im)< 0: interpbg[i] = interpbg[i-1]##if the image is all negative, use the previous value just to get past this 
               # pdb.set_trace()
                if np.nanmax(im) > 0: 
                    girdded = interpolate.griddata(okcoords2,imim[okcoords2[:,0],okcoords2[:,1]],(imyg,imxg),method=bgtype)
                    interpbg[i] = np.mean(girdded[xin,yin])

                if i == 0: 
                    bgimage = girdded*1.0/(len(hdu[1].data)*1.0)
                else: bgimage += girdded*1.0/(len(hdu[1].data)*1.0)
                    
                if i % 100 == 0: print('Interpolating BG for frame ' + str(np.max([1,i])) + ' out of ' + str(len(interpbg)))
           
            for jj in range(20):
                lowmark = np.nanpercentile(flatim[nonzero],(jj+1)*5)
                lowpix = np.where((flatim[nonzero]<lowmark) & (np.isnan(flatim[nonzero]) == False))[0]
                medim[i,jj] = np.nanmean(flatim[nonzero][lowpix])

            
            
            
            greater = np.nanpercentile(flatim[nonzero],70)
            lowpix = np.where((flatim[nonzero]<greater) & (np.isnan(flatim[nonzero]) == False))[0]
            bg[i] = np.nanmedian(flatim[nonzero[lowpix]])
            #pdb.set_trace()

        else:
            fluxcirc[i] = np.sum(im*cshape)/np.sum(cshape)
            xcen[i] = tot_image.shape[0]/2+1
            ycen[i] = tot_image.shape[1]/2+1  
            bg[i] = np.sum(im*rshape)/np.sum(rshape)
                 
        if plotmovie == True:
            plt.pcolormesh(xim,yim,im,cmap='cubehelix')
            plt.plot(xcen[i],ycen[i],'r.')
            plt.plot(xim[xin,yin],yim[xin,yin],'r.')
            pdb.set_trace()
            
            
  
    if sector == 3: ##something wrong with S3 where quals are messed up somewhat and offset by one, fix with hack here
        zxc = np.where(qual != 0)[0]
        qual[zxc-1] = qual[zxc]
        
    okok = np.where(qual==0)[0] 
    ##apply a pointing correction?
    if centroidmode == 'free':
        ##calibrate the centroid
        xbinval,binsize,ebin,allgoodindx,xbin = lcbin(cadtime,xcen,50,usemean=False,userobustmean=True,linfit=False)
        ybinval,binsize,ebin,allgoodindy,ybin = lcbin(cadtime,ycen,50,usemean=False,userobustmean=True,linfit=False)        
        usexcen = np.interp(cadtime,cadtime[allgoodindx[1:]],xcen[allgoodindx[1:]])
        useycen = np.interp(cadtime,cadtime[allgoodindy[1:]],ycen[allgoodindy[1:]])

        xoff= usexcen-np.median(usexcen)
        yoff = useycen-np.median(useycen)
        offset = (xoff**2+yoff**2)**0.5
        normoff = offset/np.median(offset[okok])
    
    ##pick the best medim BG to sub:
    ##abd regions from after datadumps
    #desat = np.where(qual == 2048)[0]
    #badmask = np.zeros(len(cadtime),dtype=int)
    #for dd in desat:
    #    bad = np.where((cadtime-cadtime[dd] < 2.0) & (cadtime-cadtime[dd] >= 0.0))[0]
    #    badmask[bad] = 1
    #okok = np.where((qual == 0) & (badmask == 0 ) & (offset < 0.05))[0]
    
    
    
    cdp = np.zeros(20)
    for jj in range(20): cdp[jj] =  np.std(fluxcirc[okok]-medim[okok,jj])
    best = np.argmin(cdp)
    bg = medim[:,best]
    if bgtype =='percentile': print('Using lowest ' + str(5*(best+1)) + ' percent of pixels as the sky level')
    if bgtype =='simple': 
        print('Using median of non-star pixels as sky level')
        bg=standard_medim*1.0
    
    ##diagnosticplot lcs, commented out ususally     
    #lc1 = (fluxcirc[okok]-interpbg[okok])/np.median(fluxcirc[okok]-interpbg[okok])
    #lc2 =(fluxcirc[okok]-bg[okok])/np.median(fluxcirc[okok]-bg[okok])
    #lc3 = (fluxcirc[okok]-standard_medim[okok])/np.median(fluxcirc[okok]-standard_medim[okok])
    #pdb.set_trace()
    
    
    if (bgtype == 'linear') | (bgtype == 'cubic'):
        bg = interpbg*1.0
##PLOTTING AND OUTPUTS
    
    if noplots == False:
        ooo = resultfilename.split('/')
        nnn = ['/'.join(ooo[0:len(ooo)-1])+'/',ooo[-1].split('.')[0]]

        plot1 = nnn[0]+'plots/'+nnn[1]+'_totimage.pdf'
        plot2 = nnn[0]+'plots/'+nnn[1]+'_maskimage.pdf'
        fig,ax = plt.subplots()
        ii = ax.pcolormesh(np.log10(tot_image-np.nanmin(tot_image)+1),cmap='cubehelix')    
        xin,yin = np.where(cshape >0)
        ax.plot(yin+0.5,xin+0.5,'r.')   
        fig.colorbar(ii,ax=ax)
        fig.savefig(plot1)
    
        fig2,ax2 = plt.subplots()
        ii2=ax2.pcolormesh(cshape,cmap='cubehelix')
        fig2.colorbar(ii2,ax=ax2)
        fig2.savefig(plot2)
    
    
        fig3,ax3 = plt.subplots()
        okok = np.where(qual==0)[0]
        ax3.plot(cadtime[okok],fluxcirc[okok]-np.nanmedian(fluxcirc),'b.',label='Raw')
        ax3.plot(cadtime[okok],(fluxcirc[okok]-bg[okok])-np.nanmedian(fluxcirc-bg),'r.',label='BG-sub')
        ax3.plot(cadtime[okok],(fluxcirc[okok]-medim[okok,0])-np.nanmedian(fluxcirc-medim[:,0]),'kx',label='BGold-sub')
        ax3.set_xlabel('Time')             
        ax3.set_ylabel('Rel. Flux')
        ax3.legend()
        plot3 = nnn[0]+'plots/'+nnn[1]+'_LCplot.pdf'
        fig3.savefig(plot3)
        
        if (bgtype =='linear') | (bgtype=='cubic'):
            fig4,ax4 = plt.subplots()
            ii=ax4.pcolormesh(np.log10(bgimage-np.nanmin(bgimage)+1),cmap='cubehelix')
            fig4.colorbar(ii,ax=ax4)
            fig4.savefig(nnn[0]+'plots/'+nnn[1]+'_bgimage.pdf')
            plt.close(fig4)    
        #plt.show()
        #pdb.set_trace()
        plt.close(fig)
        plt.close(fig2)     
        plt.close(fig3)  

    #pdb.set_trace()
    ##now make the lightcurve output structure
    output = np.recarray((len(fluxcirc),),dtype=[('t',float),('xcen',float),('ycen',float),('bg',float),('flx',float),('circrad',float),('quality',int)])
    output.t       = cadtime
    output.xcen    = xcen
    output.ycen    = ycen
    output.bg      = bg
    output.flx     = fluxcirc
    output.flx_err = fluxcirc_err
    output.quality = qual
    
    if contig == False: output.circrad = xsig*5
    #pdb.set_trace()
    if contig == True: output.circrad = -1
    ##pdb.set_trace()
    if nosave==False:
        fff = open(resultfilename,'wb')
        pickle.dump(output,fff)
        fff.close()
    if return_image==True: 
        return output,tot_image,cshape
    else: return output
    
def LCconvert(lcdata,sector3=False,correctmdump=True):
    
    ##convert to notch pipeline format inputs:
    ##is the sector is 3, there's something funky with the quality flags, make sure to fix this 
    dl = len(lcdata)
    outdata         = np.recarray((dl,),dtype=[('t',float),('fraw',float),('fcor',float),('s',float),('qual',int),('detrend',float),('divisions',float),('fraw_err',float),('fcor_err',float)])
    outdata.t = lcdata.t
    outdata.fraw = lcdata.flx
    outdata.fcor = lcdata.flx-lcdata.bg
    outdata.qual = lcdata.quality
    outdata.fraw_err = lcdata.flx_err
    outdata.fcor_err = lcdata.flx_err

    okok = np.where(np.isnan(outdata.fcor)==False)[0]  
    outdata = outdata[okok]
    outdata.fraw_err /= np.nanmedian(outdata.fraw)
    outdata.fcor_err /= np.nanmedian(outdata.fcor)
    outdata.fraw /= np.nanmedian(outdata.fraw)
    outdata.fcor /= np.nanmedian(outdata.fcor)
    
    ##remove bad points
    #pdb.set_trace()
    ##level out the momentum dump flux shifts. This is at qual==32 flag.
    ##basic idea: fit a line on either side of jump to a 0.5 day window and adjust by offset at jump time
    ###correctmdump=False
    if correctmdump==True:
        md = np.where(outdata.qual==32)[0]
        #orig=outdata.fcor*1.0
        for dd in md:
            drng=np.where((outdata.t<outdata.t[dd]+1.) & (outdata.t>outdata.t[dd]-1) & (outdata.qual == 0))[0]
            qwe1=np.where(outdata.t[drng] < outdata.t[dd])[0]
            qwe2=np.where(outdata.t[drng] > outdata.t[dd])[0]
            #pdb.set_trace()
            ###only do this if ther are enough points, and not a strong imbalance of points on the right/left of the break. This usually menas edge of a obs run and no point correcting here.
            if (len(qwe2)>10) & (len(qwe1) >10) & ((1+len(qwe1))*1.0/(1+len(qwe2)) > 0.5) & ((1+len(qwe2))*1.0/(1+len(qwe1)) > 0.5):
                lmodt,lmod,lres,qwe1,lpars,lerr=clipfitline(outdata.t[drng],outdata.fcor[drng],outdata.t[dd],niter=5,nojump=False)
                fmod = np.polyval(np.polyfit(outdata.t[drng[qwe1]],outdata.fcor[drng[qwe1]],2),outdata.t[drng[qwe1]])
                fres = (fmod-outdata.fcor[drng[qwe1]])/lerr
                lres2=(lmod[qwe1]-outdata.fcor[drng[qwe1]])/lerr
                shiftline = -lpars[3]
                bikf =np.sum(fres**2)+3*np.log(len(qwe1))                
                bik2    = np.sum(lres**2)+5*np.log(len(qwe1))
#                 print(bikf-bik2)
#                 print(shiftline)
                ##check there's really a jump based on gradient
                grad = np.gradient(outdata.fcor[drng],outdata.t[drng])
                #grad2 = np.gradient(outdata.fcor[drng])
                mg = np.median(grad)
                mad  = np.median(np.absolute(grad-mg))
                rms = np.sqrt(np.mean((grad-np.mean(grad))**2))
                #pdb.set_trace()
                #plt.plot(outdata.t[drng],outdata.fcor[drng],'xk')
               # plt.plot(outdata.t[drng[qwe1]],outdata.fcor[drng[qwe1]],'gx')
                if (np.absolute(grad[qwe2[0]-1]-mg) > mad/0.67449*3) & (bikf-bik2>20.0):outdata.fcor[dd:]+= shiftline*1.0
                #plt.plot(outdata.t[drng[qwe1]],outdata.fcor[drng[qwe1]],'.r')
               # plt.plot(lmodt,lmod,'r')
               # plt.plot(lmodt[qwe1],fmod,'g')
               # plt.show()
               # if bikf-bik2>20: pdb.set_trace()
    
    if sector3 == False: keep = np.where(outdata.qual == 0)[0]
    if sector3 == True: ##the sectorthree funky case
        badmask = np.zeros(len(outdata))
        bb = np.where(outdata.qual != 0)[0]
        bb2 = bb-1
        badmask[bb]= 1
        badmask[bb2] =1 
        keep = np.where(badmask ==0)[0]
            
    outdata = outdata[keep]
    outdata.s = 0.0
    return outdata

def LCconvertCTL(lcdata):
    ##convert to notch pipeline format inputs:
    
    dl = len(lcdata)
    outdata         = np.recarray((dl,),dtype=[('t',float),('fraw',float),('fcor',float),('s',float),('qual',int),('detrend',float),('divisions',float)])
    outdata.t = lcdata['TIME']
    outdata.fraw = lcdata['SAP_FLUX']
    outdata.fcor = lcdata['PDCSAP_FLUX']
    outdata.qual = lcdata['QUALITY']
    okok = np.where(np.isnan(outdata.fcor)==False)[0]
    outdata.fraw /= np.nanmedian(outdata.fraw)
    outdata.fcor /= np.nanmedian(outdata.fcor)
    outdata = outdata[okok]
    keep = np.where(outdata.qual == 0)[0]
    outdata.s = 0.0
    outdata = outdata[keep]
    return outdata    
    
def LCcombine(lclist):
    ##combine multiple sectors of data:
    nsect = len(lclist)
    from numpy.lib.recfunctions import stack_arrays
    outdata = stack_arrays(lclist,asrecarray= True,usemask = False)
    return outdata
    
    
    
    
    
    
