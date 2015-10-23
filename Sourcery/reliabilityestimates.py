# Reliability estimation script
# Lerato Sebokolodi <mll.sebokolodi@gmail.com>



import matplotlib
matplotlib.use('Agg')
import utils
import numpy 
import tempfile
import Tigger 
import pylab

import os
import math

class load(object):


    def __init__(self, imagename, psfname=None, sourcefinder_name='pybdsm',
                 makeplots=True, do_psf_corr=True, do_local_var=True,
                 psf_corr_region=2, local_var_region=10, rel_excl_src=None, 
                 pos_smooth=1.6, neg_smooth=1.6, loglevel=0, thresh_pix=5,
                 thresh_isl=3, neg_thresh_isl=3, neg_thresh_pix=5,
                 prefix=None, do_rel=False, **kw):

        """ Takes in image and extracts sources and makes 
            reliability estimations..
           
 
        imagename: Fits image
        psfname: PSF fits image, optional. 

        sourcefinder_name: str, optional. Default 'pybdsm'.
            Uses source finder specified by the users.

        makeplots: bool, optional. Default is True.
            Make reliability plots.

        do_psf_corr : bool, optional. Default True.
            If True, correlation of sources with PSF will be added
            as an extra source parameter in reliability estimation.
            But the PSF fits image must be provided.

        do_local_var : bool, optional. Default is True.
            Adds local variance as an extra source parameter,
            similar to do_psf_corr but independent of the PSF image. 

        psf_corr_region : int, optional. Default value is 2. 
            Data size to correlate around a source in beam sizes.
 
        local_var_region: int, optional. Default 10.
            Data size to compute the local variance in beam sizes.

        rel_excl_src : float numbers, optional. Default is None. 
            Excludes sources in this region from the reliability
            estimations, e.g ra, dec, radius in degrees. For
            many regions: ra1, dec1, radius1: ra2, dec2, radius2.

        pos_smooth : float, optional. Default 1.6
            Data smoothing threshold in the positive side of an image.
            For default value 1.6, data peaks < 1.6 * image noise
            will be averaged out.

        neg_smooth : float, optional. Default 1.6.
            Similar to pos_smooth but applied to the negative side of
            an image.

        loglevel :  int, optional. Default is 0.
            Provides Pythonlogging options, 0, 1, 2 and 3 for info, debug,
            error and critial respectively.

        thresh_isl :  float, optional. Default is 3.
            Threshold for the island boundary in number of sigma above
            the mean. Determines extent of island used for fitting 
            [pybdsm]. For positive pixels.

        thresh_pix : float, optional. Default is 5.
            Source detection threshold: threshold for the island 
            peak in number of sigma above the mean. For positive pixels.

        neg_thresh_isl : float, optional. Default is 3. 
            Simialr to thresh_isl but applied to negative side 
            of the image.

        neg_thresh_pix : float, optional. Default is 5. 
            Similar to thresh_pix but applied to the negative
            side of an image.
  
        do_rel: boolean. Default is False. 
            If True then a only reliable sources will be catalogued.
        
         kw : kward for source extractions. Should be a mapping e.g
            kw['thresh_isl'] = 2.0 or kw['do_polarization'] = True 
        """


        # log level  
        self.loglevel = loglevel
        self.log = utils.logger(self.loglevel)
       
        # image, psf image
        self.imagename = imagename
        self.psfname = psfname 
        self.log.info("Loading Image data")

        # reading imagename data
        self.imagedata, self.wcs, self.header, self.pixelsize =\
            utils.reshape_data(self.imagename)

        if not self.psfname:
            self.log.info("No psf provided, do_psf_corr = False.")
            self.do_psf_corr = False
     
        # computing negative noise
        self.noise = utils.negative_noise(self.imagedata)
        
        self.log.info("The negative noise is %e"%self.noise)

        if self.noise == 0: 
            self.log.debug("The negative noise is 0, check image")

        # source finder initialization
        self.sourcefinder_name  = sourcefinder_name
        self.log.info("Using %s source finder to extract sources."%
                      self.sourcefinder_name)

        # setting output file names       
        self.prefix = prefix
        self.poslsm = self.prefix + "_positive.lsm.html"
        self.neglsm = self.prefix + "_negative.lsm.html"


        # making negative image
        self.negativeimage = utils.invert_image(
                               self.imagename, self.imagedata,
                               self.header, self.prefix)

        # boolean optionals    
        self.makeplots = makeplots       
        self.do_psf_corr = do_psf_corr
        self.do_local_var = do_local_var
        self.do_rel = do_rel

        # smoothing factors
        self.pos_smooth = pos_smooth
        self.neg_smooth = neg_smooth
        
        # region to evaluate
        self.psf_corr_region = psf_corr_region
        self.local_var_region = local_var_region
        self.rel_excl_src = rel_excl_src
 
        # Pybdsm or source finder fitting thresholds
        self.thresh_isl = thresh_isl
        self.thresh_pix = thresh_pix
        self.opts_pos = dict(thresh_pix=self.thresh_pix,
                             thresh_isl=self.thresh_isl)
        self.opts_pos.update(kw)
        self.opts_neg = {}
        self.neg_thresh_isl = neg_thresh_isl
        self.neg_thresh_pix = neg_thresh_pix
        self.opts_neg["thresh_isl"] = self.neg_thresh_isl
        self.opts_neg["thresh_pix"] = self.neg_thresh_pix



    def source_finder(self, image=None, thresh=None,
                      noise=None, lsmname=None, **kw):
        
        #TODO look for other source finders and how they operate
         

        thresh = thresh or self.pos_smooth
        image = image or self.imagename

        ext = utils.fits_ext(image)
        tpos = tempfile.NamedTemporaryFile(suffix="."+ext, dir=".")
        tpos.flush()
        
        # data smoothing
        mask, noise = utils.thresh_mask(
                          image, tpos.name,
                          thresh=thresh, noise=self.noise, 
                          sigma=True, smooth=True)

        lsmname = lsmname or self.poslsm
        # source extraction
        utils.sources_extraction(
             image=tpos.name, output=lsmname, 
             sourcefinder_name=self.sourcefinder_name, 
             blank_limit=self.noise/100.0, **kw)


    def remove_sources_within(self, catalog, rel_excl_src=None):
   
        model = Tigger.load(catalog)
        sources = model.sources
        
        if rel_excl_src:
            for i in range(len(rel_excl_src)):
                ra, dec, tolerance = rel_excl_src[i].split(",")
                ra, dec, tolerance = map(numpy.deg2rad, (float(ra),
                                         float(dec), float(tolerance)))
                within = model.getSourcesNear(ra, dec, tolerance)   
                for src in sorted(sources):
                    if src in within:
                         sources.remove(src)
            model.save(catalog)

    
    def params(self, sources):

        labels = dict(size=(0, "Log$_{10}$(Source area)"), 
                      peak=(1, "Log$_{10}$( Peak flux [Jy] )"), 
                      tot=(2, "Log$_{10}$( Total flux [Jy] )"))
                     
        nsrc = len(sources)
        if self.do_psf_corr:
            labels.update( {"coeff":(len(labels),
                            "Log$_{10}$ (CF)")})
        if self.do_local_var:
            labels.update( {"local": (len(labels),
                            "Log$_{10}$(Local Variance)")})
        
        out = numpy.zeros([nsrc, len(labels)])

        for i, src in enumerate(sources):
            # get source extent in arcsec
            try:
               ex = numpy.rad2deg(src.get_attr("_pybdsm_Maj")) * 3600 
            except AttributeError:
                ex = bmaj * 3600
            ex = ex or bmaj * 3600
            try:
                ey = numpy.rad2deg(src.get_attr("_pybdsm_Min")) * 3600 
            except AttributeError:
                ey = bmin * 3600
            ey = ey or bmin * 3600
            area = ex * ey * math.pi

            if self.do_local_var:
                local_variance = src.l
            if self.do_psf_corr:
                cf = src.cf
            flux = src.brightness()
            peak = src.get_attr("_pybdsm_Peak_flux")
            
            if self.do_psf_corr and self.do_local_var:
                out[i,...] = area, peak, flux, cf, local_variance
            elif self.do_psf_corr:
                out[i,...] = area, peak, flux , cf
            elif self.do_local_var:
                out[i,...] = area, peak, flux , local_variance           
            else:
                out[i,...] = area, peak, flux
        return numpy.log10(out), labels 


    def get_reliability(self):


        # finding sources 
        self.source_finder(image=self.imagename, lsmname=self.poslsm, 
                           thresh=self.pos_smooth, **self.opts_pos)

        self.source_finder(image=self.negativeimage, lsmname=self.neglsm,
                           thresh=self.neg_smooth, **self.opts_neg)

        # removing sources within a specified radius
        self.remove_sources_within(catalog=self.poslsm, rel_excl_src=
                                   self.rel_excl_src)
        self.remove_sources_within(catalog=self.neglsm, rel_excl_src=
                                   self.rel_excl_src)

        # add local variance as a parameter
        if self.do_local_var:
            utils.local_variance(self.imagedata, self.header, 
                              catalog=self.poslsm, wcs=self.wcs, 
                              pixelsize=self.pixelsize, local_region=
                              self.local_var_region, savefig=False,
                              highvariance_factor=None, neg_side=True)

            utils.local_variance(self.imagedata, self.header,
                              catalog=self.neglsm, wcs=self.wcs,
                              pixelsize=self.pixelsize, local_region=
                              self.local_var_region, savefig=False,
                              highvariance_factor=None, neg_side=True)
        # compute correlation if only do_psf_corr = True 
        #and the psf is provided 
        if self.do_psf_corr and self.psfname:
            utils.psf_image_correlation(
                 catalog=self.poslsm, psfimage=self.psfname,
                 imagedata=self.imagedata, header=self.header,
                 wcs=self.wcs, pixelsize=self.pixelsize,
                 corr_region=self.psf_corr_region)
            utils.psf_image_correlation(
                 catalog=self.neglsm, psfimage=self.psfname, 
                 imagedata=self.imagedata, header=self.header,
                 wcs=self.wcs, pixelsize=self.pixelsize, 
                 corr_region=self.psf_corr_region)
      
        ##TODO verbose vs. logging
        pmodel = Tigger.load(self.poslsm, verbose=self.loglevel)
        nmodel = Tigger.load(self.neglsm, verbose=self.loglevel)
        
        posSources = pmodel.sources
        negSources = nmodel.sources

        npsrc = len(posSources)
        nnsrc = len(negSources)      
 
        positive, labels = self.params(posSources)
        negative, labels = self.params(negSources)

        # setting up a kernel, Gaussian kernel
        bandwidth = []
        for plane in positive.T:
            bandwidth.append(plane.std())

        nplanes = len(labels)
        cov = numpy.zeros([nplanes, nplanes])

        for i in range(nplanes):
            for j in range(nplanes):
                if i == j:
                    cov[i, j] = bandwidth[i]*(4.0/((nplanes+2)*
                                   npsrc))**(1.0/(nplanes+4.0))
        
        pcov = utils.gaussian_kde_set_covariance(positive.T, cov)
        ncov = utils.gaussian_kde_set_covariance(negative.T, cov)
    
        # get number densities
        nps = pcov(positive.T) * npsrc
        nns = ncov(positive.T) * nnsrc

        # define reliability of positive catalog
        rel = (nps-nns)/nps
    
        for src, rf in zip(posSources, rel):
            src.setAttribute("rel", rf)
            out_lsm = self.poslsm
        pmodel.save(out_lsm)

        if self.makeplots:
            savefig = self.prefix + "_planes.png"
            utils.plot(positive, negative, rel=rel, labels=labels,
                        savefig=savefig)

        # setting up the reliability threshold
        # get number densities
        nnps = pcov(negative.T) * npsrc
        nnns = ncov(negative.T) * nnsrc
     
        nrel = (nnps-nnns)/nnps
        reliable = nrel.max()
        self.log.info("Reliable sources have reliability > %.3f"
                      %reliable)
        if reliable < 0:
            reliable = 0.0
        if self.do_rel:
            os.system("tigger-convert --select='rel>%.3f' %s %s -f"
                      %(reliable+0.1,self.poslsm,self.poslsm))
        return  self.poslsm, self.neglsm      

