# Reliability estimation script
# Lerato Sebokolodi <mll.sebokolodi@gmail.com>



import matplotlib
matplotlib.use('Agg')
import utilss
import numpy 
import tempfile
import Tigger 
import pylab
import datetime
import os
import math

class compute(object):


    def __init__(self, imagename, psfname=None, sourcefinder_name='pybdsm',
                 makeplots=True, do_psf_corr=True, do_local_var=True,
                 psf_corr_region=2, local_var_region=10, local_negside=True, 
                 rm_sources_region=None, pos_smooth=1.6, neg_smooth=1.6,
                 loglevel=0, neg_thresh_isl=3, neg_thresh_pix=5, 
                 outdir=None, **kw):

        """ Takes in image and extracts sources and makes relibibility estimation
            for each extracted source.
           
 
            imagename: Fits image
            psfname: PSF fits image, optional. 
            sourcefinder_name: str, optional. Default 'pybdsm'.
                 Uses source finder of the users choice.
            makeplots: bool, optional. Default is True.
                 If True the reliability density plots will be made.
            do_psf_corr : bool, optional. Default True.
                 If True, correlation of sources with PSF will be added
                 as an extra source parameter in making reliability estimates.
                 But the PSF fits image must be provided.
            do_local_var : bool, optional. Default is True.
                 If True then local variance will be added as an extra source
                 parameter and will be used for making reliability estimates. 
            psf_corr_region : int, optional. Default value is 2. 
                 This specifies a factor of PSF sizes to correlate e.g.,
                 for a default value then the region size is given as
                 2 * beam size.
            local_var_region: int, optional. Default 10.
                 This specifies a factor of PSF sizes to compute local variance
                 e.g., for a default value then the region size is given as
                 10 * beam size.
            rm_sources_region : str of numbers in units of degrees, optional.
                 Default is None.
                 This gives a region in an image to which the sources inside
                 the radius will not be included for reliability estimations e.g.,
                 rm_sources_region = 'ra,dec,radius' where ra,dec and radius are
                 in degrees.
            pos_smooth : float between 0 and 1, optional. Default 0.8.
                 This smoothens the positive side on a Fits data in order to
                 make faint artefact and source peaks to be significant.
                 With large values implying more smoothing. 
            neg_smooth : float between 0 and 1, optional. Default 0.8.
                 This smoothens the negative side on a Fits data in order to
                 make faint artefact peaks to be significant over noise peaks.
            loglevel :  int, optional. Default is 0.
                 Provides logging options, 0 is info, 1 debug, 2 error and 3 
                 critial.
            neg_thresh_isl : float, optional. Default is 3. 
                 This sets source finder island threshold for the negative pixels
                 in an image.
            neg_thresh_pix : float, optional. Default is 5. 
                 This sets source finder peak threshold for the negative pixels
                 in an image.
            outdir : Will dump all output here. Default is 'sourcery-<time stamp>'
            kw : kward for source extractions. Should be a mapping e.g
                kw['thresh_isl'] = 2.0 or kw['do_polarization'] = True 
        """

        self.outdir = outdir or os.path.join("Sourcery"+datetime.datetime.now().\
                                strftime("%Y-%m-%d-%H-%M-%S"))

        self.directory = os.mkdir(self.outdir)
         
        self.loglevel = loglevel
        self.log = utilss.logger(self.loglevel)
        
        self.imagename = imagename
        self.psfname = psfname
        
        
        self.log.info("Loading Image data")
        self.imagedata, self.wcs, self.header, self.pixelsize =\
                                utilss.reshape_data(self.imagename)

        if not self.psfname:
            self.log.info("No psf provided setting do_psf_corr to False.")
            self.do_psf_corr = False
     
 
        self.noise = utilss.negative_noise(self.imagedata)
        
        self.log.info("The negative noise is %e"%self.noise)
        if self.noise == 0: 
            self.log.debug("The negative noise is 0, check image")

        self.sourcefinder_name  = sourcefinder_name
        self.log.info("Using %s source finder to extract sources."%
                     self.sourcefinder_name)

        # making negative image
        self.negativeimage = utilss.invert_image(self.imagename,
                             self.imagedata, self.header)

       
        self.do_psf_corr = do_psf_corr
        self.do_local_var = do_local_var
        self.negside = local_negside
        self.rm_sources_region = rm_sources_region

        self.pos_smooth = pos_smooth
        self.neg_smooth = neg_smooth
  
        self.psf_corr_region = psf_corr_region
        self.local_var_region = local_var_region
        self.makeplots = makeplots       
 
        
        self.opts_pos = {}
        self.opts_pos.update(kw)
        
        self.opts_neg = {}
        self.neg_thresh_isl = neg_thresh_isl
        self.neg_thresh_pix = neg_thresh_pix
        self.opts_neg["thresh_isl"] = self.neg_thresh_isl
        self.opts_neg["thresh_pix"] = self.neg_thresh_pix

    def source_finder(self, image=None, thresh=None, noise=None, **kw):
        
        #TODO look for other source finders and how they operate
         

        thresh = thresh or self.pos_smooth
        image = image or self.imagename

        ext = utilss.fits_ext(image)
        tpos = tempfile.NamedTemporaryFile(suffix="."+ext, dir=".")
        tpos.flush()
        
        mask, noise = utilss.thresh_mask(image, tpos.name,
                      thresh=thresh, noise=self.noise, sigma=True, smooth=True)
        self.out = image.replace("."+ext, ".lsm.html")

        utilss.sources_extraction(image=tpos.name, output=self.out, 
                    sourcefinder_name=self.sourcefinder_name, 
                    blank_limit=self.noise/100.0, **kw)

        return self.out
        

    def remove_sources_within(self, catalog,
                  rm_sources=None):
   
        model = Tigger.load(catalog)
        sources = model.sources
        
        if rm_sources:
            for i in range(len(rm_sources)):
                ra, dec, tolerance = rm_sources[i].split(",")
                ra, dec, tolerance = map(numpy.deg2rad, (float(ra),
                    float(dec), float(tolerance)))
                within = model.getSourcesNear(ra, dec, tolerance)   
                for src in sorted(sources):
                    if src in within:
                         sources.remove(src)
            model.save(catalog)
    
    def params(self, sources, do_coef=None, do_local_variance=None):

        labels = dict(size=(0, "Log$_{10}$(Source area)"), 
                 peak=(1, "Log$_{10}$( Peak flux [Jy] )"), 
                 tot=(2, "Log$_{10}$( Total flux [Jy] )"))
                     
        nsrc = len(sources)
        if do_coef:
            labels.update( {"coeff":(len(labels),
                  "Log$_{10}$ (CF)")})
        if do_local_variance:
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

            if do_local_variance:
                local_variance = src.l
            if do_coef:
                cf = src.cf
            flux = src.brightness()
            peak = src.get_attr("_pybdsm_Peak_flux")
            
            if do_coef and do_local_variance:
                out[i,...] = area, peak, flux, cf, local_variance
            elif do_coef:
                out[i,...] = area, peak, flux , cf
            elif do_local_variance:
                out[i,...] = area, peak, flux , local_variance           
            else:
                out[i,...] = area, peak, flux
        return numpy.log10(out), labels 
     

        
    def get_reliability(self):


        # finding sources 
        pos_catalog = self.source_finder(image=self.imagename,
                      thresh=self.pos_smooth, **self.opts_pos)
        neg_catalog = self.source_finder(image=self.negativeimage,
                       thresh=self.neg_smooth, **self.opts_neg)

        # removing sources within a specified radius
        self.remove_sources_within(catalog=pos_catalog, rm_sources=
                 self.rm_sources_region)
        self.remove_sources_within(catalog=neg_catalog, rm_sources=
                 self.rm_sources_region)

        # add local variance as a parameter
        utilss.local_variance(self.imagedata, self.header, catalog=pos_catalog,
                 wcs=self.wcs, pixelsize=self.pixelsize, local_region=
                 self.local_var_region, savefig=False, highvariance_factor=None,
                 neg_side=self.negside)
        utilss.local_variance(self.imagedata, self.header, catalog=neg_catalog,
                 wcs=self.wcs, pixelsize=self.pixelsize, local_region=
                 self.local_var_region, savefig=False, highvariance_factor=None,
                 neg_side=self.negside)
        
        if self.do_psf_corr:
            utilss.psf_image_correlation(catalog=pos_catalog, psfimage=self.psfname,
                    imagedata=self.imagedata, header=self.header, wcs=self.wcs, 
                    pixelsize=self.pixelsize, corr_region=self.psf_corr_region)
            utilss.psf_image_correlation(catalog=neg_catalog, psfimage=self.psfname, 
                    imagedata=self.imagedata, header=self.header, wcs=self.wcs,
                    pixelsize=self.pixelsize, corr_region=self.psf_corr_region)
       
        pmodel = Tigger.load(pos_catalog)
        nmodel = Tigger.load(neg_catalog)
        
        posSources = pmodel.sources
        negSources = nmodel.sources

        npsrc = len(posSources)
        nnsrc = len(negSources)      
 
        positive, labels = self.params(posSources, do_coef=self.do_psf_corr,
                          do_local_variance=self.do_local_var)
        negative, labels = self.params(negSources, do_coef=self.do_psf_corr,
                          do_local_variance=self.do_local_var)
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
        
        pcov = utilss.gaussian_kde_set_covariance(positive.T, cov)
        ncov = utilss.gaussian_kde_set_covariance(negative.T, cov)
    
        # get number densities
        nps = pcov(positive.T) * npsrc
        nns = ncov(positive.T) * nnsrc

        # define reliability of positive catalog
        rel = (nps-nns)/nps
    
        for src, rf in zip(posSources, rel):
            src.setAttribute("rel", rf)
            out_lsm = pos_catalog
        pmodel.save(out_lsm)

        if self.makeplots:
            savefig = out_lsm.replace(".lsm.html","_rel.png")
            utilss.plot(positive, negative, rel=rel, labels=labels,
                       savefig=savefig)
        self.log.info(">>>> Finding a reliable value of the reliability <<<<")
        # setting up the reliability threshold
        # get number densities
        nnps = pcov(negative.T) * npsrc
        nnns = ncov(negative.T) * nnsrc
     
        nrel = (nnps-nnns)/nnps
        reliable = nrel.max()
        self.log.info("Reliable sources have reliability > %.3f"%reliable)

