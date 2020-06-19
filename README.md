# THYMEffi
Code to make FFI light curves, with background correction and dynamically sized apertures from TESS full frame images for the TESS Hunt for Young and Maturing Exoplanets (THYME) team and others.

Setup and usage is described in the scripts, but feel free to ask Aaron Rizzuto for help at any time.

# What this code actually does:
1. Downloads a tesscut TPF from the tesscut website api.
2. Take the median image accross the campaign, and fit a Gaussian to determine the aperture size to use for extraction, there are ways of forcing this to be a certain values if you want to wth fixaperture=X (pixels). It's also limited to not be larger than 7 pixels or less than 3 pixels.   
3. Then goes frame by frame doing aperture photometry, with radius of target aperture as above. Computes the background using one of three methods (bgtype='simple' or 'percentile' or 
linear'  or 'cubic') see below
4. Corrects background and formats the output to something mostly standardized for ACR's codes.


Background methods:
set with e.g. bgtype = 'percentile' in call to run_extraction

1. 'percentile'
Chooses for the whole lightcurve which percentile to use for the sky level based on standard deviation across lightcurve after subtraction of each background level. For HIP 67522, that turned out to be the 40% lowest pixels being used as the sky pixels. Taking the median of the lowest x percentile of pixels as the sky level for percentiles of 5-100 in steps of 5%.

2. 'simple' 
    Just takes the median pixel value across the image, this only really works for very bright things.

3. 'linear' and 'cubic': 
    Automatically masks most stars from the image, the target included, and then interpolates either a linear or cubic function over the remaining image and masked pixels. 
    This seems to give the best results for the fainter things, or cases where percentile leaves systematics. Still not perfect though.!
    Note: this has trouble with the triple-stripe feature in sector 11 (and elsewhere?), ACR still working on this. Also has trouble when there's too much crowding and not enough to interpolate from.
    Note that when one of these modes is set, the plot directory produced by run_extraction for each target extracted will contain a plot called *******_bgimage.pdf. Which will show the mean background interpolated image 
    for the sector. Worth checking to make sure it did't bork on something and interpolate crazily. 
    

Adapted from Aaron Rizzuto's development repository, will be updated as important changes happen.
Does not include the Notch filter pipeline, though the lightcurves produced from this are compatible with those codes (see https://bitbucket.org/aaronrizzuto/detrend-k2 for that)
