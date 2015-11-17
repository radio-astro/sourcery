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
                 prefix=None, do_nearsources=False, increase_beam_cluster=False,
                 **kw):

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

        do_nearsources: boolean. Default is False.
            If true it adds number of nearest neighnours as an extra
            parameter. It looks for sources around 5 beam sizes.

        savefits: boolean. Default is False.
            if True a negative image is saved.

        increase_beam_cluster: boolean, optional. If True, sources
            Gaussian groupings will be increase by 20%. If False,
            the actual beam size will be used. Default is False.
   
         kw : kward for source extractions. Should be a mapping e.g
            kw['thresh_isl'] = 2.0 or kw['do_polarization'] = True 
        """


       
        # image, psf image
        self.imagename = imagename
        self.psfname = psfname 
        # setting output file names  
     
        self.prefix = prefix
        self.poslsm = self.prefix + ".lsm.html"
        self.neglsm = self.prefix + "_negative.lsm.html"

        # log level  
        self.loglevel = loglevel
        self.log = utils.logger(self.loglevel, prefix=self.prefix)

        self.log.info(" Loading the image data")

 

        # reading imagename data
        self.imagedata, self.wcs, self.header, self.pixelsize =\
            utils.reshape_data(self.imagename, prefix=self.prefix)

        self.bmaj = numpy.deg2rad(self.header["BMAJ"])

        self.do_psf_corr = do_psf_corr
        if not self.psfname:
            self.log.info(" No psf provided, do_psf_corr is set to False.")
            self.do_psf_corr = False
           
        # computing negative noise
        self.noise = utils.negative_noise(self.imagedata)
        
        self.log.info(" The negative noise is %e Jy/beam"%self.noise)

        if self.noise == 0: 
            self.log.debug(" The negative noise is 0, check image")

        # source finder initialization
        self.sourcefinder_name  = sourcefinder_name
        self.log.info(" Using %s source finder to extract the sources."%
                      self.sourcefinder_name)


        # making negative image
        self.savefits = False
        self.negativeimage = utils.invert_image(
                               self.imagename, self.imagedata,
                               self.header, self.prefix)

        
        # increase the beam by 20% on its major and minor axis
        self.do_beam = increase_beam_cluster
        
        # boolean optionals    
        self.makeplots = makeplots
        self.do_local_var = do_local_var
        self.nearsources = do_nearsources

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
        self.bmin, self.bpa =  self.header["BMIN"], self.header["BPA"]

        if self.do_beam:
            bmaj = self.header["BMAJ"]
            self.opts_pos["beam"] = (1.2*bmaj, 1.2*self.bmin, self.bpa)

        
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
             blank_limit=self.noise/100.0, prefix=self.prefix,
             **kw)


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

    def nearest_neighbour(self, model, src):
        near = model.getSourcesNear(src.pos.ra, src.pos.dec, 5 * self.bmaj)
        return (1.0/float(len(near)))
    
    def params(self, sources, model=None):

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
        if self.nearsources:
            labels.update( {"Nearer": (len(labels),
                            "Log$_{10}$(Source Near)")})     
        
        out = numpy.zeros([nsrc, len(labels)])

        for i, src in enumerate(sources):
            # get source extent in arcsec
            try:
               ex = numpy.rad2deg(src.get_attr("_pybdsm_Maj")) * 3600 
            except AttributeError:
                ex = self.bmaj * 3600
            ex = ex or self.bmaj * 3600
            try:
                ey = numpy.rad2deg(src.get_attr("_pybdsm_Min")) * 3600 
            except AttributeError:
                ey = self.bmin * 3600
            ey = ey or self.bmin * 3600
            area = ex * ey * math.pi

            if self.do_local_var:
                local_variance = src.l
            if self.do_psf_corr:
                cf = src.cf

            if self.nearsources:
                near = self.nearest_neighbour(model, src)

            flux = src.brightness()
            peak = src.get_attr("_pybdsm_Peak_flux")
  
            if self.do_psf_corr and self.do_local_var and self.nearsources:
                out[i,...] = area, peak, flux, cf, local_variance, near

            elif self.do_psf_corr and self.nearsources: 
                out[i,...] = area, peak, flux, cf, near

            elif self.do_local_var and self.nearsources: 
                out[i,...] = area, peak, flux, local_variance, near

            elif self.do_psf_corr and self.do_local_var:
                out[i,...] = area, peak, flux, cf, local_variance

            elif self.do_psf_corr:
                out[i,...] = area, peak, flux , cf

            elif self.do_local_var:
                out[i,...] = area, peak, flux , local_variance 

            elif self.nearsources:
                out[i,...] = area, peak, flux , near                     
            else:
                out[i,...] = area, peak, flux

        return numpy.log10(out), labels 


    def get_reliability(self):


        # finding sources 
        self.log.info(" Extracting the sources on both sides ")
        
        self.source_finder(image=self.imagename, lsmname=self.poslsm, 
                           thresh=self.pos_smooth, **self.opts_pos)

        self.source_finder(image=self.negativeimage, lsmname=self.neglsm,
                           thresh=self.neg_smooth, **self.opts_neg)

        self.log.info(" Source Finder completed successfully ")

        if not self.savefits:
            os.system("rm -r %s"%self.negativeimage)

        # removing sources within a specified radius
        
        self.remove_sources_within(catalog=self.poslsm, rel_excl_src=
                                   self.rel_excl_src)
        self.remove_sources_within(catalog=self.neglsm, rel_excl_src=
                                   self.rel_excl_src)

        # add local variance as a parameter
        if self.do_local_var:
            self.log.info(" Computing the local variance around positive detections ")
            utils.local_variance(self.imagedata, self.header, 
                              catalog=self.poslsm, wcs=self.wcs, 
                              pixelsize=self.pixelsize, local_region=
                              self.local_var_region, savefig=False,
                              highvariance_factor=None, prefix=self.prefix,
                              neg_side=True)

            self.log.info(" DONE: Local variance on the positive side was sucessful ")
            self.log.info(" Computing the local variance around the negative detections ")

            utils.local_variance(self.imagedata, self.header,
                              catalog=self.neglsm, wcs=self.wcs,
                              pixelsize=self.pixelsize, local_region=
                              self.local_var_region, savefig=False,
                              highvariance_factor=None, prefix=self.prefix, neg_side=True)

            self.log.info("DONE: Computation of the local variance completed successfully ")
           
        # compute correlation if only do_psf_corr = True 
        #and the psf is provided 
        if self.do_psf_corr and self.psfname:
            self.log.info(" Computing the correlation factor of the detections with the PSF ")
            utils.psf_image_correlation(
                 catalog=self.poslsm, psfimage=self.psfname,
                 imagedata=self.imagedata, header=self.header,
                 wcs=self.wcs, pixelsize=self.pixelsize,
                 corr_region=self.psf_corr_region, prefix= self.prefix)
            utils.psf_image_correlation(
                 catalog=self.neglsm, psfimage=self.psfname, 
                 imagedata=self.imagedata, header=self.header,
                 wcs=self.wcs, pixelsize=self.pixelsize, 
                 corr_region=self.psf_corr_region, prefix=self.prefix)
            self.log.info(" DONE: Correlation factor has been successfully assigned to the detections.")
        ##TODO verbose vs. logging
        pmodel = Tigger.load(self.poslsm, verbose=self.loglevel)

        nmodel = Tigger.load(self.neglsm, verbose=self.loglevel)
        
        posSources = pmodel.sources
        negSources = nmodel.sources

        npsrc = len(posSources)
        nnsrc = len(negSources)      
 
        positive, labels = self.params(posSources, pmodel)
        negative, labels = self.params(negSources, nmodel)

        # setting up a kernel, Gaussian kernel
        bandwidth = []

        for plane in negative.T:
            bandwidth.append(plane.std())



        nplanes = len(labels)
        cov = numpy.zeros([nplanes, nplanes])


        for i in range(nplanes):
            for j in range(nplanes):
                if i == j:
                    cov[i, j] = bandwidth[i]*((4.0/((nplanes+2)*
                                  nnsrc))**(1.0/(nplanes+4.0)))

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
                        savefig=savefig, prefix=self.prefix)

        return  self.poslsm, self.neglsm      

