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


def logger(level=0):
    logging.basicConfig()

    LOGL = {"0": "INFO",
            "1": "DEBUG",
            "2": "ERROR",
            "3": "CRITICAL"}

    log = logging.getLogger(" Sourcery ")
    log.setLevel(eval("logging."+LOGL[str(level)]))

    return log
#-----------------------knicked from Stats.py------------------------------- 
def reshape_data (image):

    """ Reshape FITS data to (stokes,freq,npix_ra,npix_dec).
        Returns reshaped data, wcs, the image header,
        and pixel size 
         
        image: Fits data  
    """

    with pyfits.open(image) as hdu:
        data = hdu[0].data
        hdr = hdu[0].header
        shape = list(data.shape)
        ndim = len(shape)

    wcs = WCS(hdr, mode='pyfits')
    
    pixel_size = abs(hdr['CDELT1'])

    if ndim<2:
        raise ValueError('The FITS file needs at least two dimensions')

 # This is the shape I want the data in
    want = (
            ['STOKES',0],
            ['FREQ',1],
            ['RA',2],
            ['DEC',3],
)
   
    # Assume RA,DEC is first (FITS) or last two (NUMPY)
    if ndim>3:
        for ctype,ind in want[:2]:
            for axis in range(1,ndim+1):
                if hdr['CTYPE%d'%axis].startswith(ctype):
                    want[ind].append(ndim-axis)
        if want[0][-1] == want[1][-2] and want[0][-2] == want[1][-1]:
            tmp = shape[0]
            shape[0] = shape[1]
            shape[1] = tmp
            data = numpy.reshape(data,shape)
    if ndim ==3:
        if not hdr['CTYPE3'].startswith('FREQ'):
            data = data[0,...]
    elif ndim>4:
        raise ValueError('FITS file has more than 4 axes. Aborting')
        
    return data,wcs,hdr,pixel_size


def negative_noise(data):

    """ Computes the image noise using the negative pixels """

    negative = data[data<0].flatten()
    noise = numpy.concatenate([negative,-negative]).std()
    return noise

def invert_image(image, data, header):
    
    ext = fits_ext(image)
    output = image.replace(ext,'_negative.fits')
    newdata = -data
    pyfits.writeto(output, newdata, header, clobber=True)
    return output


def thresh_mask(imagename, outname, thresh, noise=None, sigma=False, smooth=None):
    """ Create a threshhold mask """

    hdu = pyfits.open(imagename)
    hdr = hdu[0].header
    
    ndim = hdr["NAXIS"] 
    imslice = [0]*ndim
    imslice[-2:] = [slice(None)]*2
    data = hdu[0].data[imslice].copy()

    # If smooth is not specified, use a fraction of the beam
    
    if sigma:
        noise = noise or negative_noise(data)
    else:
        noise = 1

    thresh = thresh*noise
    
    mask = numpy.ones(data.shape)

    if smooth:
        emin = hdr["BMIN"]
        emaj = hdr["BMAJ"]
        cell = abs(hdr["CDELT1"])
        beam = math.sqrt(emin*emaj)/cell
        scales = [.1, .33, .5, 1., 1.5, 2.]#, 4., 8., 16.]
        smooth = None
        for scale in scales: 
            kk = scale*beam
            smooth = filters.gaussian_filter(data if smooth is None else smooth, [kk,kk])
            mask *= smooth<thresh
    else:
        mask = data<thresh
    
    hdu[0].data *= (mask==False)
    hdu.writeto(outname,clobber=True)

    return mask==False, noise


def sources_extraction(image, output=None, sourcefinder_name='pybdsm', **kw):


    """Runs pybdsm on the specified 'image', converts the 
       results into a Tigger model and writes it to output.

    image :  Fits image data
    output : Tigger format, default image name.lsm.html
           A Catalog name to store the extracted sources
    """
    
    ext = fits_ext(image)
    output = output or image.replace(ext, '.lsm.html')
    gaul = output+'.gaul'
    # start with default PYBDSM options
    opts = {}
    opts.update(kw)
     
    if sourcefinder_name.lower() == 'pybdsm':
        from lofar import bdsm
        img = bdsm.process_image(image,group_by_isl=True,**kw)
        img.write_catalog(outfile=gaul,format='ascii',catalog_type='gaul',clobber=True)
    verifyGaulModel(gaul)

    # converting the model to Tigger
    tc = ['tigger-convert',gaul,output,"-t","Gaul","-f","--rename",'-o','Tigger']

    process = subprocess.Popen([' '.join(['%s'%item for item in tc])],
                  stderr=subprocess.PIPE if not isinstance(sys.stderr,
                        file) else sys.stderr,
                  stdout=subprocess.PIPE if not isinstance(sys.stdout,
                        file) else sys.stdout,
                  shell=True)

    if process.stdout or process.stderr:
        out,err = process.comunicate()
        sys.stdout.write(out)
        sys.stderr.write(err)
        out = None;
    else:
        process.wait()
    if process.returncode:
        print 'tigger-convert returns errr code %d'%(process.returncode)
    else:
        print 'DONE: tigger-convert succeeded. catalog is at %s'%output



#---------------------------------------------------------------------------------------
#knicked from pyxis lsm.pybdsm_search
def verifyGaulModel(gaullsm):
  """Check all sources in a gaul file are in valid locations before running tigger
  convert. Useful when images are 'all-sky' and have undefined regions.
  """
  falseSources = 0
  olsm = ''
  names = []
  fh=open(gaullsm, 'r')
  for ll in fh.readlines():
    cll = ' '.join(ll.split())
    if cll == '' or cll.startswith('#'):
      olsm += ll
      if cll.startswith("# Gaus_id"):
        names = cll.split()
      continue
    lineArray = cll.split(' ')
    if math.isnan(float(lineArray[names.index("RA")] )) : falseSources += 1
    if float(lineArray[names.index("Peak_flux")]) <= 0 : 
        falseSources+=1
    if float(lineArray[names.index("Total_flux")]) <= 0 :
        falseSources += 1
    else: olsm += ll
  fh.close()

  fh=open(gaullsm, 'w')
  fh.write(olsm)
  fh.close() 
     


#------------------------------------------------------- 
class gaussian_kde_set_covariance(stats.gaussian_kde):
    def __init__(self, dataset, covariance):
        self.covariance = covariance
        stats.gaussian_kde.__init__(self, dataset)
    def _compute_covariance(self):
        #if numpy.linalg.det(self.covariance) != 0:
        self.inv_cov = pylab.linalg.inv(self.covariance)
        self._norm_factor = numpy.sqrt(pylab.linalg.det(2*numpy.pi*self.covariance)) * self.n
#---------------------------------------------------------------------------


def fits_ext(fitsname):
    ext = fitsname.split(".")[-1]
    return ext


def local_variance(imagedata, header, catalog, wcs, pixelsize, tag=None,
                   local_region=5, noise=None, savefig=True, 
                   highvariance_factor=0.8, localvariance_tag=None,
                   neg_side=False):

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
        local_region: int, optional. A default value of 5. 
              Gives a region to compute the local variance in
              psf sizes, e.g, 'local_region = 2',
              then a region (= 2 * beam size) around a source is used.
        highvariance_factor: float, optional. A default value of 0.8. 
              If highvariance_factor=0.8 is given this means that
              the sources with local variance greater than  
              0.8*image_noise will be tagged 'high_variance' if
              localvariance_tag=None.
        localvariance_tag : str, optional. A default tag 'high_variance'
              A tag provided to sources of high variance as
              determined by high_variance factor. 
    """

    data = imagedata
    beam = header['BMAJ']

    bmaj = int(round(beam/pixelsize)) # beam size in pixels
    
    step = local_region * bmaj
    noise = noise or negative_noise(data)
    
    model = Tigger.load(catalog, verbose=False)
    sources = []

    if tag: 
        sources = filter(lambda src: src.getTag(tag), model.sources) 
    else:
        for src in model.sources:
            sources.append(src)
    
    positions_sky = [map(lambda rad: numpy.rad2deg(rad),\
                    (src.pos.ra,src.pos.dec)) for src in sources]
    positions = [wcs.wcs2pix(*pos) for pos in positions_sky]

    shape = data.shape 
    ndim = len( data.shape)
    if ndim == 4:
        data = data[0,0,...]
    if ndim == 3:
        data = data[0,...]
 
    step = [step,step]
    
     
    for i, (pos, src) in enumerate(zip(positions, sources)):
        x,y = pos
        if x>shape[-2] or y>shape[-1] or numpy.array(pos).any()<0:
            positions.remove(pos)
            model.sources.remove(src)
            sources.remove(src)

        if (y+step[1] > shape[-1]) or (y-step[1] < 0):
            if pos in positions:
                positions.remove(pos)
                model.sources.remove(src)
                sources.remove(src)

        if (x+step[0] > shape[-2]) or (x-step[0] < 0):
            if pos in positions:
                positions.remove(pos)
                model.sources.remove(src)
                sources.remove(src)
    _std = []
    
    if neg_side:
        data = -data
 
    n = 0
    for (x, y), srs in zip(positions,sources):
        pos = [x,y]
        subrgn = data[y-step[0]:y+step[0], x-step[1]:x+step[1]]
        subrgn = subrgn[subrgn > 0]
        std = subrgn.std()
        if math.isnan(float(std)) or _std == 0:
            sources.remove(srs)
            model.sources.remove(srs)
            positions.remove(srs)
            n += 1
        else:
            _std.append(std)
            srs.setAttribute('l', std)
        
    if n > 0:
        print("Nan encountered %d times. Increase the size of the\
               region or check the image. Otherwise sources with 0 or nan\
               are flagged." %n)


    def high_variance_sources(pos, local_variance, noise, model, threshold,
                savefig=savefig, localtag=None):

        if savefig:
            save_fig = catalog.replace('.lsm.html','.png')

        x = numpy.arange(len(pos))
        pylab.figure()
        pylab.plot(x, local_variance)
        pylab.plot([noise] * len(local_variance))
        localtag = localtag or 'high_variance'    
        for i,(pos,src) in enumerate(zip( pos, model.sources)):
            if _std[i] > threshold:
                src.setTag(localtag, True)
                pylab.plot(x[i], local_variance[i], 'rD')
                pylab.annotate(src.name, xy=(x[i],local_variance[i]))
        if savefig:
            pylab.ylabel('local variance')
            pylab.savefig(save_fig)
    

    if highvariance_factor:
        threshold = highvariance_factor * noise
        high_variance_sources(positions, _std, noise, model, threshold=threshold,
                 savefig=savefig, localtag=localvariance_tag)
    model.save(catalog)   
    return _std 


def psf_image_correlation(catalog, psfimage, imagedata, header, wcs ,
                     pixelsize, corr_region=5, thresh=0.4, tags=None,
                     coefftag=None, do_high=False):


    """ Computes correlation between the image and PSF image

        catalog : Source model, Tigger format.
        psfimage : Instrument's Point spread functions Fits data.
        imagedata : Fits data
        header : Fits header e.g., img=pyfits.open("test.fits")
              header=img[0].header
        wcs : Obtained from Fit header e.g.,
              wcs = WCS(header,mode='pyfits') where WCS is imported as
              'from astLib.astWCS import WCS'.
              This class provides methods for accessing information from the World
              Coordinate System (WCS) contained in the header of a FITS image.
              Conversions between pixel and WCS coordinates can also be performed.
        pixelsize: float, obtained in degrees from Fits header e.g 
              header['BMAJ'], where BMAJ is beam major axis.
        corr_region : int, optional. A default value of 5.
              corr_region sets a region in Fits data and psf data to be correlated, e.g.,
              for a default value this implies that a region = 5* beam size 
              will be correlated. 
        thresh : float, optional. A default value of 0.4. 
              Gives a correlation threshold i.e those sources
              with correlation > thresh will be tagged as high
              correlation sources using a tag 'high_corr' if coefftag=None.
        tags: str, optional. 
              If tag is provided then sources which has the specified tag will be
              correlated with the psf.
        coefftag: str, optional. A Default string is 'high_corr'.
              If provided sources with correlation > thresh will be tagged
              using the user specified tag. 
        do_high: bool, optional.  
              If True, sources of high correlation are tagged using 'coefftag',
              if False no tagging will be made.
    """

    model = Tigger.load(catalog,verbose=False)
   
    image_data = imagedata 
    beam = header['BMAJ']
    psf_data, wcs_psf, psf_hdr, psf_pix = reshape_data(image=psfimage)
    
    shape = image_data.shape
 
    bmaj = int(round(beam/pixelsize))
    step = corr_region * bmaj
    
    sources = []
    if tags: 
        sources = filter(lambda src: src.getTag(tags),model.sources) 
    else:
         for src in model.sources:
             sources.append(src)
    
    positions_sky = [map(lambda rad: numpy.rad2deg(rad),
                    (src.pos.ra,src.pos.dec))  for src in sources]
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

    for i,(p,src) in enumerate(zip(pos,sources)):      
        x,y = p
        if x>shape[-2] or y>shape[-1] or numpy.array(p).any()<0:
            pos.remove(p)
            sources.remove(src)
            model.sources.remove(src)

        if (y+step[1] > shape[-1]) or (y-step[1] < 0):
            if p in pos:
                pos.remove(p)
                model.sources.remove(src)
                sources.remove(src)
        if (x+step[0] > shape[-2]) or (x-step[0] < 0):
            if p in pos:
                pos.remove(p)
                model.sources.remove(src)
                sources.remove(src)

    central = psf_hdr['CRPIX2']
    psf_region = psf_data[central-step[0] : central+step[0],
                 central-step[1] : central+step[1]]
    psf_region = psf_region.flatten()
     
    corr = []
    n = 0
    for src,(ra,dec) in zip(sources,pos): 
        data_region = image_data[dec-step[0] : dec+step[0],
                      ra-step[1] : ra+step[1]].flatten()
        norm_data = (data_region-data_region.min())/(data_region.max()-
                    data_region.min())
        c_region = numpy.corrcoef((norm_data,psf_region))
        cf_region =  (numpy.diag((numpy.rot90(c_region))**2).sum())**0.5/2**0.5
        cf = cf_region

        if math.isnan(float(cf)) or cf == 0:
            model.sources.remove(src)
            sources.remove(src)
            n += 1
        else:
            src.setAttribute('cf',cf)
            corr.append(cf)

    if n > 0:
        print('%d sources were removed due to 0 or nan correlation'%n)

    thresh = thresh 
    coefftag = coefftag or 'high_corr'
    if do_high:
        for src,crr in zip(sources,corr):
            if crr > thresh:
               src.setTag(coefftag,True)
    model.save(catalog)     
    return corr

def plot(pos, neg, rel=None, labels=None, show=False, savefig=None):
    #info('Making Reliability plots')
 
    # labels for projections
    plots = []
    nplanes = len(labels)
    for i, label_i in enumerate(labels):
        for label_j in labels.keys()[i+1:]:
            i = labels[label_i][0]
            j = labels[label_j][0]
            plots.append( [i, j, label_i, label_j] )

    
    npos, nneg = len(pos), len(neg)
    pylab.figure(figsize=(10*nplanes, 8*nplanes))

    if nneg < 5:
        print 'Warn here' ##TODO        
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

        pylab.subplot( int(row), int(column), counter+1)
        a,b = neg[:, i], neg[:, j]
        c,d = pos[:, i], pos[:, j]

        kernel = numpy.array([a.std(),b.std()])
        cov = numpy.array([(kernel[0]**2,0.0),(0.0,kernel[1]**2)])
        ncov = gaussian_kde_set_covariance(numpy.array([a,b]),cov)
        
        # define axis limits for plots
        ac = numpy.concatenate((a,c))        
        bd = numpy.concatenate((b,d))        
        pylab.xlim(ac.min(), ac.max()*1.2)
        pylab.ylim(bd.min(), bd.max()*1.2)

        #negative detection density field
        PN = ncov(numpy.array([a,b])) * nneg
       
        xi = numpy.linspace(ac.min() ,ac.max(), 100)
        yi = numpy.linspace(bd.min(), bd.max(), 100)
        zzz = griddata((a,b),PN,(xi[None,:],yi[:,None]),method='cubic')

        pylab.contour(xi,yi,zzz,20,linewidths=4) 
        pylab.scatter(pos[:,i],pos[:,j],marker='o',c=rel,s=40)
        pylab.xlabel(labels[x][1],fontsize='35')
        pylab.ylabel(labels[y][1],fontsize='35')
        pylab.grid()

    if savefig : 
        pylab.savefig(savefig)
    if show : 
        pylab.show()
    pylab.clf()

