# THYMEffi
Code to make FFI light curves, with background correction and dynamically sized apertures from TESS full frame images for the TESS Hunt for Young and Maturing Exoplanets (THYME) team and others.

Setup and usage is described in the scripts, but feel free to ask Aaron Rizzuto for help at any time.

# What this code actually does:
1. Download a tesscut TPF from the tesscut website api.
2. Take the median image accross the campaign, and fit a Gaussian to determine the aperture size to use for extraction, there are ways of forcing this to be a certain values if you want to wth fixaperture=X (pixels). It's also limited to not be larger than 7 pixels or less than 3 pixels.  
3. Then goes frame by frame doing aperture photometry, with radius of target aperture as above, taking the medium of the lowest x percentile of pixels as the sky level for percentiles of 5-100 in steps of 5%.
4. Chooses for the whole lightcurve which percentile to use for the sky level based on standard deviation across lightcurve after subtraction of each background level. For HIP 67522, that turned out to be the 40% lowest pixels being used as the sky pixels. 

Adapted from Aaron Rizzuto's development repository, will be updated as important changes happen.
Does not include the Notch filter pipeline, though the lightcurves produced from this are compatible with those codes (see https://bitbucket.org/aaronrizzuto/detrend-k2 for that)
