# Reliability estimator and direction-dependent source tagging functions tools



import matplotlib
matplotlib.use('Agg')
import Tigger
import pyfits
import numpy
import subprocess
import tempfile
import os 
import sys
import logging
import pylab
from scipy.ndimage import filters
from astLib.astWCS import WCS
import math
from scipy import stats
from scipy.interpolate import griddata


matplotlib.rcParams.update({'font.size': 12})

def logger(level=0):
    
    name = "logfile.log"
    logging.basicConfig(filename=name)

    LOGL = {"0": "INFO",
            "1": "DEBUG",
            "2": "ERROR",
            "3": "CRITICAL"}

    log = logging.getLogger(" Sourcery ")
    log.setLevel(eval("logging."+LOGL[str(level)]))

    return log


log = logger(level=0)


#-----------------------knicked from Stats.py------------------------------- 
def reshape_data (image):

    """ Reshape FITS data to (stokes,freq,npix_ra,npix_dec).

    Returns reshaped data, wcs, the image header, and
    pixel size.
         
    image: Fits data  
    """

    with pyfits.open(image) as hdu:
        data = hdu[0].data
        hdr = hdu[0].header
        shape = list(data.shape)
        ndim = len(shape)

    wcs = WCS(hdr, mode="pyfits")
    
    pixel_size = abs(hdr["CDELT1"])

    if ndim < 2:
        log.error("The FITS file needs at least two dimensions")
        

 # This is the shape I want the data in
    want = (
            ["STOKES",0],
            ["FREQ",1],
            ["RA",2],
            ["DEC",3],
)
   
    # Assume RA,DEC is first (FITS) or last two (NUMPY)
    if ndim > 3:
        for ctype, ind in want[:2]:
            for axis in range(1, ndim+1):
                if hdr["CTYPE%d"%axis].startswith(ctype):
                    want[ind].append(ndim-axis)
        if want[0][-1] == want[1][-2] and want[0][-2] == want[1][-1]:
            tmp = shape[0]
            shape[0] = shape[1]
            shape[1] = tmp
            data = numpy.reshape(data,shape)
    if ndim == 3:
        if not hdr["CTYPE3"].startswith("FREQ"):
            data = data[0,...]
    elif ndim > 4:
        log.error("FITS file has more than 4 axes. Aborting")
        
    return data, wcs, hdr, pixel_size


def negative_noise(data):

    """ Computes the image noise using the negative pixels """

    negative = data[data<0].flatten()
    noise = numpy.concatenate([negative,-negative]).std()

    return noise


def invert_image(image, data, header, prefix=None):
    
    ext = fits_ext(image)
    output = prefix + "_negatives.fits" or image.replace(ext,"_negative.fits")
    newdata = -data
    pyfits.writeto(output, newdata, header, clobber=True)

    return output


def thresh_mask(imagename, outname, thresh, 
                noise=None, sigma=False, smooth=None):
    """ Create a threshhold mask """

    hdu = pyfits.open(imagename)
    hdr = hdu[0].header
    
    ndim = hdr["NAXIS"] 
    imslice = [0] * ndim
    imslice[-2:] = [slice(None)] * 2
    data = hdu[0].data[imslice].copy()

    # If smooth is not specified, use a fraction of the beam
    
    if sigma:
        noise = noise or negative_noise(data)
    else:
        noise = 1

    thresh = thresh * noise
    
    mask = numpy.ones(data.shape)

    if smooth:
        emin = hdr["BMIN"]
        emaj = hdr["BMAJ"]
        cell = abs(hdr["CDELT1"])
        beam = math.sqrt(emin*emaj)/cell
        scales = [.1, 2, 5.0, 10., 20, 40.]#, 4., 8., 16.]
        smooth = None
        for scale in scales: 
            kk = scale * beam
            smooth = filters.gaussian_filter(
                       data if smooth is None else 
                       smooth, [kk,kk])
            mask *= smooth < thresh
    else:
        mask = data < thresh
    
    hdu[0].data *= (mask==False)
    hdu.writeto(outname, clobber=True)

    return mask==False, noise


def sources_extraction(image, output=None,
                       sourcefinder_name="pybdsm",
                       **kw):


    """Runs pybdsm on the specified 'image', converts the 
       results into a Tigger model and writes it to output.

    image :  Fits image data
    output : Tigger format, default image name.lsm.html
        A Catalog name to store the extracted sources
    """
    
    ext = fits_ext(image)
    output = output or image.replace(ext, ".lsm.html")
    gaul = output + ".gaul"
    # start with default PYBDSM options
    opts = {}
    opts.update(kw)
     
    if sourcefinder_name.lower() == "pybdsm":
        from lofar import bdsm
        img = bdsm.process_image(image, group_by_isl=True, **kw)
        img.write_catalog(outfile=gaul, format="ascii", 
                          catalog_type="gaul", clobber=True)
    verifyGaulModel(gaul)

    # converting the model to Tigger
    tc = ["tigger-convert", gaul, output,
          "-t","Gaul","-f","--rename","-o","Tigger"]

    process = subprocess.Popen([" ".join(["%s"%item for item in tc])],
                  stderr = subprocess.PIPE if not isinstance(sys.stderr,
                        file) else sys.stderr,
                  stdout = subprocess.PIPE if not isinstance(sys.stdout,
                        file) else sys.stdout,
                  shell=True)

    if process.stdout or process.stderr:
        out,err = process.comunicate()
        sys.stdout.write(out)
        sys.stderr.write(err)
        out = None
    else:
        process.wait()
    if process.returncode:
        log.error("tigger-convert returns errr code %d"%
                 (process.returncode))
    else:
        log.info("DONE: tigger-convert succeeded catalog is %s"%output)



#---------------------------------------------------------------------------------------
#knicked from pyxis lsm.pybdsm_search
def verifyGaulModel(gaullsm):

  """Check all sources in a gaul file are in valid locations
     before running tigger.

  convert. Useful when images are 'all-sky' and have undefined regions.
  """

  falseSources = 0
  olsm = ""
  names = []
  fh=open(gaullsm, "r")
  for ll in fh.readlines():
    cll = ' '.join(ll.split())
    if cll == '' or cll.startswith("#"):
      olsm += ll
      if cll.startswith("# Gaus_id"):
        names = cll.split()
      continue
    lineArray = cll.split(" ")
    if math.isnan(float(lineArray[names.index("RA")] )) : 
        falseSources += 1
    if float(lineArray[names.index("Peak_flux")]) <= 0 : 
        falseSources+=1
    if float(lineArray[names.index("Total_flux")]) <= 0 :
        falseSources += 1
    else: olsm += ll
  fh.close()

  fh=open(gaullsm, "w")
  fh.write(olsm)
  fh.close() 
     


#----------------------------------------------------
#knicked from Sofia reliability estimator
class gaussian_kde_set_covariance(stats.gaussian_kde):
    def __init__(self, dataset, covariance):
        self.covariance = covariance
        stats.gaussian_kde.__init__(self, dataset)
    def _compute_covariance(self):
        #if numpy.linalg.det(self.covariance) != 0:
        self.inv_cov = pylab.linalg.inv(self.covariance)
        self._norm_factor = numpy.sqrt(
             pylab.linalg.det(2*numpy.pi*self.covariance)) * self.n
#---------------------------------------------------------------------------


def fits_ext(fitsname):
    ext = fitsname.split(".")[-1]
    return ext


def local_variance(imagedata, header, catalog, wcs=None,
                   pixelsize=None, tag=None,local_region=5,
                   noise=None, savefig=True, highvariance_factor=0.8,
                   high_local_tag=None, neg_side=False, 
                   setatr=True, do_high_loc=False, prefix=None):

    """ Calculates the local varience (lv) around a source on 
        one side of interest. 
 
    imagedata : Reshaped Fits data
    header : Image header
    catalog : Source catalog, in Tigger format.
        Source model to compute local variance around.

    tag : str, optional.
        if specified then the local variance will be
        computed for only a subset of sources with a tag,
        e.g., 'tag=snr'.
    header : Fits header e.g., img=pyfits.open("test.fits")
        header=img[0].header
    wcs :This class provides methods
        for accessing information from the World  Coordinate System
        (WCS) contained in the header of a FITS image. Conversions 
        between pixel and WCS coordinates can also be performed.
        If not provided it is directly obtained from the Fits header
        provided.
    pixelsize: float, Default is None.
         If not provided then it is directly obtained form a Fits header.
    local_region: int, optional. A default value of 5. 
        Gives a region to compute the local variance in
        psf sizes, e.g, 'local_region = 2',
        then a region (= 2 * beam size) around a source is used.

    highvariance_factor: float, optional. A default value of 0.8. 
        If highvariance_factor=0.8 is given this means that
        the sources with local variance greater than  
        0.8*image_noise will be tagged 'high_var' if
        high_local_tag=None.

    high_local_tag : str, optional. A default tag None.
        Tag assigned to sources of high local variance.
        If None is provided the default tag will be used
        i.e 'high_var'. Else it uses the specified tag.

    setatr : bool, optional. Default is True.
        If True all sources will be tagged with 'l'
        giving each detection an extra local variance
        parameter.

    do_high_loc : bool, optional. Default is False.
        If True sources with high local variance will be tagged
        using 'localvariance_tag' (see above).
    prefix :  str, optional. Default is None.
    """

    data = imagedata
    beam = header["BMAJ"]

    if pixelsize is None:
        pixelsize = abs(header["CDELT1"])
    if wcs is None:
        wcs = WCS(header, mode="pyfits")

    bmaj = int(round(beam/pixelsize)) # beam size in pixels

    if not isinstance(local_region, int):
        if isinstance(local_region, float):
            local_region = int(round(local_region))
            log.debug("Float is provided and int is required,"
                       "arounding offto the nearest integer")
            if local_region == 0:
                log.error("It rounded off to zero now,"
                          "change local_region into an integer."
                          "Aborting")
        else:
            log.error("local_region must be an integer. Abort")
    
    step = local_region * bmaj

    noise = noise or negative_noise(data)
    
    model = Tigger.load(catalog)
    sources = []

    if tag:
        log.debug("Local variance is only computed"
                  "for sources with a tag %s are"%
                   tag)
        sources = filter(lambda src: src.getTag(tag), model.sources) 
    else:
        for src in model.sources:
            sources.append(src)
    
    positions_sky = [map(lambda rad: numpy.rad2deg(rad),
                    (src.pos.ra,src.pos.dec)) for src in sources]
    positions = [wcs.wcs2pix(*pos) for pos in positions_sky]

    shape = data.shape 
    ndim = len( data.shape)
    if ndim == 4:
        data = data[0,0,...]
    if ndim == 3:
        data = data[0,...]
 
    step = [step, step]
    
    m = 0 
    for i, (pos, src) in enumerate(zip(positions, sources)):
        x,y = pos
        if x>shape[-2] or y>shape[-1] or numpy.array(pos).any()<0:
            positions.remove(pos)
            model.sources.remove(src)
            sources.remove(src)
            m += 1

        if (y+step[1] > shape[-1]) or (y-step[1] < 0):
            if pos in positions:
                positions.remove(pos)
                model.sources.remove(src)
                sources.remove(src)
            m += 1

        if (x+step[0] > shape[-2]) or (x-step[0] < 0):
            if pos in positions:
                positions.remove(pos)
                model.sources.remove(src)
                sources.remove(src)
            m += 1
    if m > 0:
        log.debug("It will be useful to increase the image size,"
                  "sources with ra+step or dec+step > image size"
                  "are removed")
        
    _std = []
    
    if neg_side:
        data = -data
        log.debug("Using the negative side of the provided image.")
 
    n = 0
    for (x, y), srs in zip(positions, sources):
        pos = [x,y]
        subrgn = data[y-step[0] : y+step[0], x-step[1] : x+step[1]]
        subrgn = subrgn[subrgn > 0]
        std = subrgn.std()

        if math.isnan(float(std)) or _std == 0:
            sources.remove(srs)
            model.sources.remove(srs)
            positions.remove(srs)
            n += 1
        else:
            _std.append(std)
            if setatr:
                srs.setAttribute("l", std)
        
    if n > 0:
        log.debug("Nan encountered %d times. Increase the size of the"
                  "region or check the image. Otherwise sources with"
                  "0 or nan are flagged." %n)


    def high_variance_sources(
            pos, local_variance, noise, model, threshold,
            savefig=savefig, prefix=None, localtag=None):

        if savefig:
            save_fig = prefix + "_variance.png" or catalog.replace(
                                                      ".lsm.html", ".png")

        local = [l/1.0e-6 for l in local_variance]
        x = numpy.arange(len(pos))
        pylab.figure()
        pylab.plot(x, local)
        pylab.plot([noise/1.0e-6] * len(local_variance))

        localtag = localtag     
        for i, (pos, src) in enumerate(zip( pos, model.sources)):
            if _std[i] > threshold:
                src.setTag(localtag, True)
                pylab.plot(x[i], local[i], "rD")
                pylab.annotate(src.name, xy=(x[i],local[i]))
        if savefig:
            pylab.ylabel("local variance[$\mu$]")
            pylab.savefig(save_fig)

    if high_local_tag is None:
         high_local_tag = "high_var"
    if do_high_loc:
        threshold = highvariance_factor * noise
        high_variance_sources(positions, _std, noise, model,
                              threshold=threshold, savefig=savefig,
                              prefix=prefix, localtag=high_local_tag)
    model.save(catalog)   

    return _std 


def psf_image_correlation(catalog, psfimage, imagedata, header, wcs=None ,
                     pixelsize=None, corr_region=5, thresh=0.4, tags=None,
                     coefftag='high_corr', setatr=True, do_high=False):


    """ Computes correlation of the image and PSF image

    catalog : Source model, Tigger format.
    psfimage : Instrument's Point spread functions Fits data.
    imagedata : Fits data
    header : Fits header e.g., img=pyfits.open("test.fits")
        header=img[0].header
    wcs :This class provides methods
        for accessing information from the World  Coordinate System
        (WCS) contained in the header of a FITS image. Conversions 
        between pixel and WCS coordinates can also be performed.
        If not provided it is directly obtained from the Fits header
        provided.
    pixelsize: float, Default is None.
         If not provided then it is directly obtained form a Fits header.
    corr_region : int, optional. A default value of 5.
        The size of the region to correlate given in beam sizes.
    thresh : float, optional. A default value of 0.4. 
        Correlation threshold. Sources with correlation > threshold
        are sources with high correlation.
    tags: str, optional. Default is None.
        If specified only sources with 'Tag' will be evaluated.
    coefftag: str, optional. A Default string is 'high_corr'.
        If provided sources with correlation > thresh will be tagged
        using the user specified tag.
    setatr : bool, optional. Default is True.
        If True all sources will be tagged with 'cf'
        giving each detection an extra correlation with PSF parameter.
    do_high: bool, optional.  Default is False.
        If True, sources of high correlation are tagged using 'coefftag',
        if False no tagging will be made.
    """

    model = Tigger.load(catalog)
   
    image_data = imagedata 
    beam = header["BMAJ"]
    psf_data, wcs_psf, psf_hdr, psf_pix = reshape_data(image=psfimage)
    
    shape = image_data.shape

    if pixelsize is None:
        pixelsize = abs(header["CDELT1"])
    if wcs is None:
        wcs = WCS(header, mode="pyfits")

    bmaj = int(round(beam/pixelsize))
    log.info("Beam major= %d degrees"%bmaj)

    if bmaj == 0:
        log.debug("Beam major axis was read as 0, setting it to 1")
        bmaj = 1.0

    if not isinstance(corr_region, int):
        if isinstance(corr_region, float):
            corr_region = int(round(corr_region))
            log.debug("Float is provided and int is required,"
                      "arounding off to the nearest integer")
            if corr_region == 0:
                log.error("Rounding off to 0. Provide an integer."
                          "Aborting")
        else:
            log.error("corr_region must be an integer. Abort")
    
    step = corr_region * bmaj
    
    sources = []

    if tags:
        log.debug("Only sources with tag %s will be correlated"
                  "with the PSF"%tags)
        sources = filter(lambda src: src.getTag(tags),model.sources) 
    else:
         for src in model.sources:
             sources.append(src)
    
    positions_sky = [map(lambda rad: numpy.rad2deg(rad),
                    (src.pos.ra, src.pos.dec))  for src in sources]
    pos = [wcs.wcs2pix(*pos) for pos in positions_sky]

    step = [step,step]
    ndim = len(shape)
    if ndim == 4:
        image_data = image_data[0,0,...]
    if ndim == 3:
        image_data = image_data[0,...]

    pdim = len(psf_data.shape)
    if pdim == 4:
        psf_data = psf_data[0,0,...]
    if pdim == 3:
        psf_data = psf_data[0,...]
  
    
    m = 0
    for i, (p, src) in enumerate(zip(pos, sources)):      
        x,y = p
        if x>shape[-2] or y>shape[-1] or numpy.array(p).any()<0:
            pos.remove(p)
            sources.remove(src)
            model.sources.remove(src)
            m += 1

        if (y+step[1] > shape[-1]) or (y-step[1] < 0):
            if p in pos:
                pos.remove(p)
                model.sources.remove(src)
                sources.remove(src)
            m += 1
        if (x+step[0] > shape[-2]) or (x-step[0] < 0):
            if p in pos:
                pos.remove(p)
                model.sources.remove(src)
                sources.remove(src)
            m += 1
    if m > 0:
        log.debug("It will be useful to increase the image size,"
                  "sources with ra+step or dec+step > image size"
                  "are removed")

    central = psf_hdr["CRPIX2"]
    psf_region = psf_data[central-step[0] : central+step[0],
                 central-step[1] : central+step[1]]
    psf_region = psf_region.flatten()
     
    corr = []
    n = 0
    for src, (ra, dec) in zip(sources, pos): 
        data_region = image_data[dec-step[0] : dec+step[0],
                                ra-step[1] : ra+step[1]].flatten()
        norm_data = (data_region-data_region.min())/(data_region.max()-
                                                     data_region.min())
        c_region = numpy.corrcoef((norm_data, psf_region))
        cf_region =  (numpy.diag((numpy.rot90(c_region))**2)
                                  .sum())**0.5/2**0.5
        cf = cf_region

        if math.isnan(float(cf)) or cf == 0:
            model.sources.remove(src)
            sources.remove(src)
            n += 1
        else:
            corr.append(cf)

            if setatr:
                src.setAttribute("cf",cf)

    if n > 0:
        log.debug("%d sources were removed due to 0/nan correlation"%n)

    thresh = thresh 
    coefftag = coefftag
    if do_high:
        for src, crr in zip(sources, corr):
            if crr > thresh:
               src.setTag(coefftag, True)
    model.save(catalog)     

    return corr


def plot(pos, neg, rel=None, labels=None, show=False, savefig=None):

    log.info("Making Reliability plots")
 
    # labels for projections
    plots = []
    nplanes = len(labels)
    for i, label_i in enumerate(labels):
        for label_j in labels.keys()[i+1:]:
            i = labels[label_i][0]
            j = labels[label_j][0]
            plots.append( [i, j, label_i, label_j] )

    
    npos, nneg = len(pos), len(neg)
    pylab.figure(figsize=(14*nplanes, 8*nplanes))

    if nneg < 5:
        log.error("Few number of detections cant proceed plotting."
                  "Aborting")        
        return 

    if nplanes %2.0 == 0:
        row  = nplanes/2.0
        column = (nplanes-1.0)
    else:
        row =  (nplanes - 1.0)/2.0
        column = nplanes
    row_fix = row
    if  column > 4.0:
        row = column
        column = row_fix

    for counter, (i, j, x, y) in enumerate(plots):

        pylab.subplot(int(row), int(column), counter+1)
        a,b = neg[:, i], neg[:, j]
        c,d = pos[:, i], pos[:, j]

        kernel = numpy.array([a.std(), b.std()])
        cov = numpy.array([(kernel[0]**2, 0.0),(0.0, kernel[1]**2)])
        ncov = gaussian_kde_set_covariance(numpy.array([a, b]), cov)
        
        # define axis limits for plots
        ac = numpy.concatenate((a, c))        
        bd = numpy.concatenate((b, d))        
        pylab.xlim(ac.min(), ac.max()*1.2)
        pylab.ylim(bd.min(), bd.max()*1.2)

        #negative detection density field
        PN = ncov(numpy.array([a, b])) * nneg
       
        xi = numpy.linspace(ac.min() ,ac.max(), 100)
        yi = numpy.linspace(bd.min(), bd.max(), 100)
        zzz = griddata((a, b), PN,(xi[None,:], yi[:,None]), method="cubic")
        pylab.tick_params(axis='x', labelsize=30)
        pylab.tick_params(axis='y', labelsize=30)
        pylab.contour(xi, yi, zzz, 20, linewidths=4) 
        pylab.scatter(pos[:,i], pos[:,j], marker="o", c=rel, s=40)
        pylab.xlabel(labels[x][1], fontsize="35")
        pylab.ylabel(labels[y][1], fontsize="35")
        pylab.grid()

    if savefig : 
        log.info("Saving the reliability plot.")
        pylab.savefig(savefig)
    

