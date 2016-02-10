# Determing sources that require direction dependent calibration solutions
# Lerato Sebokolodi <mll.sebokolodi@gmail.com>


import matplotlib
matplotlib.use('Agg')
import Tigger
from Tigger.Coordinates import angular_dist_pos_angle as dist
import utils
import numpy
import math

class load(object):


    def __init__(self, imagename, pmodel, nmodel, header, psfname=None, noise=None,
                 snr_thresh=40, local_thresh=0.4, high_corr_thresh=0.5, negdetec_region=10, 
                 negatives_thresh=5, phasecenter_excl_radius=None,
                 prefix=None, loglevel=0):


        """ Determines sources that require direction-dependent (DD)
            calibration solutions.

        psfname: PSF fits data

        pmodel: Model of the positive image.

        nmodel: Model of the negative image
 
        header: The header of the input image

        noise: float, Default None.
             The noise of the image.

        snr_thresh: float, optional. Default is 40.
             Any source with 40 x the minimum SNR.

        local_thresh: float, optional. Default is 0.4.
             Sources with local variance > 0.4 * the noise have
             high local variance.

        high_corr_thresh:  float, optional. Default is 0.5.
             Sources of high PSF correlation have correlation above 0.5.

        negdetec_region:  float, optional. Default is 10.
             Region to lookup for negative detections around. In beam size.

        negative_thresh:  float, optional. Default is 5.
             The number of nearby negative detections. Sources
             with number > 5 require direction
             dependent (DD) calibration solutions.

        phasecenter_excl_region:  float (in degrees), optional.
             A radius from the phase center (in beam sizes) to exclude the sources
             from the evaluation.

        prefix: str, optional. Sets a prefix to the output directory.

        loglevel: int, optional. Default 0. Python logging.
                  0, 1, 2, 3 for info, debug, error and 
                  critical, respectively.
        """

        self.loglevel = loglevel
        self.prefix = prefix
        self.log = utils.logger(self.loglevel, prefix=self.prefix)
        #image and psf image

        self.pmodel = pmodel
        self.nmodel = nmodel
        self.psfname =  psfname
        self.hdr = header
  
        self.log.info(" Loading image data")

        self.noise = noise
        if not self.noise:
            self.log.info(" No noise value provided."
                          " Setting it to 1e-6. Please provide"
                          " the noise.")
            self.noise = 1e-6

        # tags
        self.dd_tag = "dE"

        # thresholds
        self.snr_thresh = snr_thresh
        self.localthresh = local_thresh
        self.corrthresh = high_corr_thresh
        self.negthresh = negatives_thresh
       
        
        #regions
        self.phaserad = phasecenter_excl_radius # radius from the phase center
        self.negregion =  negdetec_region # region to look for negatives
        
       # conversion
        self.r2d = 180.0/math.pi
        self.d2r = math.pi/180.0
        self.bmaj = self.hdr["BMAJ"] # in degrees
        
        self.ra0 =  self.hdr["CRVAL1"] * self.d2r
        self.dec0 = self.hdr["CRVAL2"] * self.d2r
 

    def number_negatives(self, source):
        

        #sources = filter(lambda src: src.getTag(tag), psources)

        tolerance = self.negregion * self.bmaj * self.d2r

        if self.phaserad:
            radius = numpy.deg2rad(self.phaserad * self.bmaj)
            
        
        ra, dec = source.pos.ra, source.pos.dec # in radians
        within = self.nmodel.getSourcesNear(ra, dec, tolerance)    
        if len(within) >= self.negthresh:
            if self.phaserad:
                if dist(self.ra0, self.dec0, ra, dec)[0] > radius: 
                        source.setTag(self.dd_tag, True)
            else:
                source.setTag(self.dd_tag, True)


    def source_selection(self):
        
        sources = self.pmodel.sources
            
        snr = [src.flux.I/self.noise for src in sources]
            
        thresh = self.snr_thresh * min(snr)
        n = 0
        for srs, s in zip(sources, snr):
            if s > thresh:
                local = srs.l 
                if local > self.localthresh * self.noise:
                    if self.psfname:
                        corr = srs.cf
                        if corr > self.corrthresh:
                            self.number_negatives(srs)

                    if not self.psfname:
                        self.number_negatives(srs)
                            
        return self.pmodel, self.nmodel
