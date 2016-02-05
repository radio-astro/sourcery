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
from Tigger.Models import SkyModel, ModelClasses
import pyfits
from Tigger.Coordinates import angular_dist_pos_angle as dist
import sys


class load(object):


    def __init__(self, imagename, psfname=None, sourcefinder_name='pybdsm',
                 makeplots=True, do_psf_corr=True, do_local_var=True,
                 psf_corr_region=5, local_var_region=10, rel_excl_src=None, 
                 pos_smooth=2, neg_smooth=2, loglevel=0, thresh_pix=5,
                 thresh_isl=3, neg_thresh_isl=3, neg_thresh_pix=5,
                 prefix=None, do_nearsources=False, increase_beam_cluster=False,
                 savemask_pos=False, savemask_neg=False, **kw):

        """ Takes in image and extracts sources and makes 
            reliability estimations..
           
 
        imagename: Fits image
        psfname: PSF fits image, optional. 

        sourcefinder_name: str, optional. Default 'pybdsm'.
            Uses source finder specified.

        makeplots: bool, optional. Default is True.
            Make reliability plots.

        do_psf_corr : bool, optional. Default True.
            If True, PSF correlation will be added
            as an extra parameter for density estimations.
            NB: the PSF fits image must be provided.

        do_local_var : bool, optional. Default is True.
            If True, adds local variance as an extra parameter,
            for density estimations. 
        
        do_nearsources: boolean. Default is False.
            If true it adds number of nearest neighnours as an extra
            parameter. It looks for sources around 5 beam sizes.

        psf_corr_region : int, optional. Default value is 5. 
            Data size to correlate around a source, in beam sizes.
 
        local_var_region: int, optional. Default 10.
            Data size to compute the local variance in beam sizes.

        rel_excl_src : floats, optional. Default is None. 
            Excludes sources in a specified region
            e.g ra, dec, radius in degrees. For
            2 regions: ra1, dec1, radius1: ra2, dec2, radius2, etc.

        pos_smooth : float, optional. Default 2.
            Masking threshold for the positive image.
            For default value 2, data peaks < 2 * image noise
            are masked.

        neg_smooth : float, optional. Default 2.
            Similar to pos_smooth but applied to the negative image.

        thresh_isl :  float, optional. Default is 3.
            Threshold for forming islands in the positive image

        thresh_pix : float, optional. Default is 5.
            Threshold for model fitting, in positive image.

        neg_thresh_isl : float, optional. Default is 3. 
            Simialr to thresh_isl but for negative image.

        neg_thresh_pix : float, optional. Default is 5. 
            Similar to thresh_pix but for negative image.

        savefits: boolean. Default is False.
            If True a negative image is saved.

        increase_beam_cluster: boolean, optional. If True, sources
            groupings will be increase by 20% the beam size. If False,
            the actual beam size will be used. Default is False.

        savemask_pos: boolean, optional. If true the mask applied on 
            the positive side of an image after smoothing is saved.
            
        savemask_neg: Similar to savemask_pos but for the negative
            side of an image.
        
        loglevel : int, optional. Default is 0.
            Provides Pythonlogging options, 0, 1, 2 and 3 are for info, debug,
            error and critial respectively.
   
         kw : kward for source extractions. Should be a mapping e.g
            kw['thresh_isl'] = 2.0 or kw['do_polarization'] = True 
        """


       
        self.prefix = prefix

        # log level  
        self.loglevel = loglevel
        self.log = utils.logger(self.loglevel, prefix=self.prefix)

        
        # image, psf image
        self.imagename = imagename
        self.psfname = psfname 
      
        # reading imagename data
        imagedata, self.wcs, self.header, self.pixelsize =\
            utils.reshape_data(self.imagename, prefix=self.prefix)

        self.imagedata = numpy.array(imagedata, dtype=numpy.float32)
        self.posdata =  numpy.array(utils.image_data(self.imagedata, self.prefix), dtype=numpy.float32)
        
        self.bmaj = numpy.deg2rad(self.header["BMAJ"])

        # boolean optionals    
        self.makeplots = makeplots
        self.do_local_var = do_local_var
        self.nearsources = do_nearsources
        self.do_psf_corr = do_psf_corr
        self.savemaskpos = savemask_pos
        self.savemaskneg = savemask_neg
        self.savefits = False

        if not self.psfname:
            self.log.info(" No psf provided, do_psf_corr is set to False.")
            self.do_psf_corr = False
        if self.psfname:
            data, self.psfhdr = utils.open_psf_image(self.psfname)
            self.psfdata = numpy.array(utils.image_data(data, self.prefix), dtype=numpy.float32)

 
        # computing negative noise
        self.noise = utils.negative_noise(self.posdata) #here is 2X2 data here
        
        self.log.info(" The negative noise is %e Jy/beam"%self.noise)
        if self.noise == 0: 
            self.log.debug(" The negative noise is 0, check image")

        # source finder initialization
        self.sourcefinder_name  = sourcefinder_name
        self.log.info(" Using %s source finder to extract the sources."%
                      self.sourcefinder_name)

        self.negimage = self.prefix + "_negative.fits"
        
        negativedata =  utils.invert_image(
                               self.imagename, self.imagedata,
                               self.header, self.negimage)

        self.negdata = numpy.array(utils.image_data(negativedata, self.prefix), numpy.float32)
        self.negativedata = numpy.array(negativedata, numpy.float32)
       
        # smoothing factors
        self.pos_smooth = pos_smooth
        self.neg_smooth = neg_smooth
        
        # region to evaluate
        self.corrstep = psf_corr_region
        self.localstep = local_var_region
        self.radiusrm = rel_excl_src
        self.do_beam = increase_beam_cluster
         

        # conversion
        self.d2r = math.pi/180.0
        self.r2d = 180.0/math.pi
        beam_pix = int(round(self.bmaj * self.r2d/self.pixelsize))
        self.locstep = self.localstep * beam_pix
        self.cfstep = self.corrstep * beam_pix
        self.bmin, self.bpa =  self.header["BMIN"], self.header["BPA"]
        if self.do_beam:
            bmaj = self.header["BMAJ"]
            self.opts_pos["beam"] = (1.2*bmaj, 1.2*self.bmin, self.bpa)

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
 

    def source_finder(self, image=None, imagedata=None, thresh=None, prefix=None,
                      noise=None, output=None, savemask=None, **kw):
        
        ext = utils.fits_ext(image)
        tpos = tempfile.NamedTemporaryFile(suffix="."+ext, dir=".")
        tpos.flush()
        
        # data smoothing
        mask, noise = utils.thresh_mask(image, 
                          imagedata, tpos.name, hdr=self.header,
                          thresh=thresh, noise=self.noise, 
                          sigma=True, smooth=True, prefix=prefix, 
                          savemask=savemask)
        naxis = self.header["NAXIS1"] 
        boundary = numpy.array([self.locstep, self.cfstep])
        trim_box = (boundary.max(), naxis - boundary.max(),
                  boundary.max(), naxis- - boundary.max())
             
        # source extraction
        utils.sources_extraction(
             image=tpos.name, output=output, 
             sourcefinder_name=self.sourcefinder_name, 
             blank_limit=self.noise/1000.0, trim_box=trim_box, prefix=self.prefix,
             **kw)
        tpos.close()


    def remove_sources_within(self, model):
        
        sources = model.sources
        for i in range(len(self.radiusrm)):
                ra, dec, tolerance = self.radiusrm[i].split(",")
                ra_r =  float(ra) * self.d2r
                dec_r =  float(dec) * self.d2r
                tolerance_r = float(tolerance) * self.d2r
                within = model.getSourcesNear(ra_r, dec_r, tolerance_r)   
                for src in sorted(sources):
                    if src in within:
                         sources.remove(src)

    
    def params(self, modelfits):
     
        # reads in source finder output             
        data = pyfits.open(modelfits)[1].data
        tfile = tempfile.NamedTemporaryFile(suffix=".txt")
        tfile.flush()
        
        # writes a catalogue in a temporaty txt file
        
        with open(tfile.name, "w") as std:
            std.write("#format:name ra_d dec_d i emaj_r emin_r pa_d\n")

        model = Tigger.load(tfile.name) # open a tmp. file

        peak, total, area = [], [], []

        for i, src in enumerate(data):
            ra, dec = data["RA"][i], data["DEC"][i]
            flux = data["Total_flux"][i] 
            emaj, emin = data["Maj"][i], data["Min"][i]
            pa = data["PA"][i]
            name = "SRC%d"%i
        
            posrd =  ModelClasses.Position(ra*self.d2r, dec*self.d2r)
            flux_I = ModelClasses.Polarization(flux, 0, 0, 0)
            shape = ModelClasses.Gaussian(emaj*self.d2r, emin*self.d2r, pa*self.d2r)
            srs = SkyModel.Source(name, posrd, flux_I, shape=shape)
           
            # area: find ex and ey if are 0 assign beam size
            if emaj or emin == 0:
                srcarea = math.pi * (self.bmaj * self.r2d) * pow(3600.0,-2) *\
                       (self.bmin *  self.r2d) 
            if  emaj and emin > 0: 
                srcarea = emaj * emin * math.pi * pow(3600.0,-2)
            

            peak_flux = data["Peak_flux"][i]

            # only accepts sources with flux > 0 and non-nan RA and DEC 
            if flux > 0 and peak_flux > 0 and not math.isnan(float(ra))\
                and not math.isnan(float(dec)):
                  model.sources.append(srs) 
                  peak.append(peak_flux)
                  total.append(flux)
                  area.append(srcarea)

        # removes sources in a given radius from the phase center
        if self.radiusrm:
            self.log.info(" Remove sources ra, dec, radius of  %r" 
                          " from the phase center" %self.radiusrm)
            self.remove_sources_within(model)


        labels = dict(size=(0, "Log$_{10}$(Source area)"), 
                      peak=(1, "Log$_{10}$( Peak flux [Jy] )"), 
                      tot=(2, "Log$_{10}$( Total flux [Jy] )"))
        if self.do_psf_corr:
            labels.update( {"coeff":(len(labels),
                            "Log$_{10}$ (CF)")})
        if self.do_local_var:
            labels.update( {"local": (len(labels),
                            "Log$_{10}$(Local Variance)")})
        if self.nearsources:
            labels.update( {"near": (len(labels),
                            "Log$_{10}$(Near Sources)")})

        nsrc = len(model.sources)
        out = numpy.zeros([nsrc, len(labels)])         
        
        # returning parameters
        for i, src in enumerate(model.sources):

            ra, dec = src.pos.ra, src.pos.dec
            pos = [self.wcs.wcs2pix(*(ra * self.r2d, dec * self.r2d))][0] #deg to pixel

            # computes the local variance
            
            loc = utils.compute_local_variance(
                        self.negdata, pos, self.locstep) 
            near = model.getSourcesNear(ra, dec, 5 * self.bmaj)
            nonear = len(near) 

            if math.isnan(float(loc)) or loc <=0:
                model.sources.remove(src)
                self.log.info(" Sources %s is removed as it has"
                              " nan or 0 local variance."%src.name)

            else:
                src.setAttribute("l", loc)
                src.setAttribute("n", nonear)
                if self.psfname:
                    correlation = utils.compute_psf_correlation(
                                  self.posdata, self.psfdata, 
                                  self.psfhdr, pos,  self.cfstep)
                    if not math.isnan(float(correlation)) or correlation >= 0:
                        corr = correlation
                        src.setAttribute("cf", corr)
                    else:
                        model.sources.remove(src)

                if self.do_psf_corr and self.do_local_var and self.nearsources:
                     out[i,...] =  area[i], peak[i], total[i], corr, loc, nonear

                elif self.do_psf_corr and self.do_local_var and not self.nearsources:
                    out[i,...] =   area[i], peak[i], total[i] , corr, loc
        
                elif self.do_psf_corr and self.nearsources and not self.do_local_var:
                    out[i,...] =   area[i], peak[i], total[i] , corr, nonear
            
                elif not self.do_psf_corr and self.do_local_var and self.nearsources:
                    out[i,...] =   area[i], peak[i], total[i] , loc, nonear
            
                elif self.do_psf_corr and not self.do_local_var and not self.nearsources:
                    out[i,...] =   area[i], peak[i], total[i] , corr
            
                elif not self.do_psf_corr and self.do_local_var and not self.nearsources:
                    out[i,...] =   area[i], peak[i], total[i] , loc
            
                elif not self.do_psf_corr and not self.do_local_var and self.nearsources:
                    out[i,...] =   area[i], peak[i], total[i] , nonear

                else:
                    out[i,...] =   area[i], peak[i], total[i]
        # removes the rows with 0s
        removezeros = (out == 0).sum(1)
        output = out[removezeros <= 0, :]
        
        return model, numpy.log10(output), labels 



    def get_reliability(self):


        # finding sources 
        self.log.info(" Extracting the sources on both sides ")
        pfile = self.prefix + ".gaul.fits"
        nfile = self.prefix + "_negative.gaul.fits"
         
        # i need to catch mmap.mmap error here

        # running a source finder
        self.source_finder(image=self.negimage, imagedata=self.negativedata,
                           output=nfile, thresh=self.neg_smooth,
                           savemask=self.savemaskneg,
                           prefix=self.prefix, **self.opts_neg)
        self.source_finder(image=self.imagename, imagedata=self.imagedata,
                           output=pfile, thresh=self.pos_smooth, 
                           savemask=self.savemaskpos,
                           prefix=self.prefix, **self.opts_pos)

        self.log.info(" Source Finder completed successfully ")
        if not self.savefits:
            self.log.info("Deleting the negative image.")
            os.system("rm -r %s"%self.negimage)

         
        pmodel, positive, labels = self.params(pfile)
        nmodel, negative, labels = self.params(nfile)
     
        # setting up a kernel, Gaussian kernel
        bandwidth = []

        for plane in negative.T:
            bandwidth.append(plane.std())

        nplanes = len(labels)
        cov = numpy.zeros([nplanes, nplanes])
        nnsrc = len(negative)
        npsrc = len(positive)

        if nnsrc == 0 or npsrc ==0:
            self.log.info("The resulting array has length of 0 thus cannot compute"
                    " the reliability. Aborting.")
        self.log.info(" There are %s positive detections "%npsrc)
        self.log.info(" There are %s negative detections "%nnsrc)

        self.log.info(" Computing the reliabilities ")
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
        for src, rf in zip(pmodel.sources, rel):
            src.setAttribute("rel", rf)
        self.log.info(" Saved the reliabilities values.")
        
        if self.makeplots:
            savefig = self.prefix + "_planes.png"
            utils.plot(positive, negative, rel=rel, labels=labels,
                        savefig=savefig, prefix=self.prefix)

        return  pmodel, nmodel, self.noise

