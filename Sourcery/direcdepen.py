# Determing sources that require direction dependent calibration solutions
# Lerato Sebokolodi <mll.sebokolodi@gmail.com>


import matplotlib
matplotlib.use('Agg')
import Tigger
from Tigger.Coordinates import angular_dist_pos_angle as dist
import utils
import numpy
import math
from astLib.astWCS import WCS

class load(object):


    def __init__(self, imagedata, psfname, pmodel, nmodel, header, local_step=10,
                 snr_thresh=40, high_corr_thresh=0.5, negdetec_region=10,
                 negatives_thresh=5, phasecenter_excl_radius=None,
                 prefix=None, loglevel=0):


        """ Determines sources that require direction-dependent (DD)
            calibration solutions.

        psfname: PSF fits data

        pmodel: Model of the positive image.

        nmodel: Model of the negative image
 
        header: The header of the input image

        snr_thresh: float, optional. Default is 40.
             Any source with 40 x the minimum SNR.

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

        # tags
        self.dd_tag = "dE"

        # thresholds
        self.snr_factor = snr_thresh
        #self.localthresh = local_thresh
        self.corrthresh = high_corr_thresh
        self.negthresh = negatives_thresh
        self.wcs = WCS(self.hdr, mode="pyfits")
        
        self.data = imagedata
        self.locstep = local_step

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
                

    def local_noise(self, pos):
        
        # computing the local noise using MAD
        x, y = pos
        
        sub_data =  self.data[y-self.locstep:y+self.locstep,
                              x-self.locstep:x+self.locstep]
        noise = numpy.mean(abs(sub_data - numpy.mean(sub_data)))
        return noise
     

    def source_selection(self):
        
        sources = self.pmodel.sources
        noise, mean = utils.negative_noise(self.data, self.prefix)
        for srs in sources:
            pos = map(lambda rad: numpy.rad2deg(rad),(srs.pos.ra,srs.pos.dec))
            positions = self.wcs.wcs2pix(*pos) # from degs to pixels 
            local_noise = self.local_noise(positions) # Local noise
            signal_to_noise = srs.flux.I/local_noise # source SNR
            thresh = self.snr_factor * noise # uses the Global Noise
            if signal_to_noise > thresh and srs.rel > 0.99:
                if self.psfname:
                    corr = srs.cf
                    if corr > self.corrthresh:
                        self.number_negatives(srs)

                if not self.psfname:
                    self.number_negatives(srs)
                            
        return self.pmodel, self.nmodel
