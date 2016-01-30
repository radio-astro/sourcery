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


    def __init__(self, imagename, poscatalog, negcatalog, psfname=None,
                 snr_thresh=100, local_thresh=0.6, local_region=10,
                 psfcorr_region=2, high_corr_thresh=0.5, negdetec_region=10, 
                 negatives_thresh=10, phasecenter_excl_radius=None,
                 prefix=None, loglevel=0):


        """ Determines sources that require direction-dependent (DD)
            calibration solutions.

        imagename: Fits data
        psfname : PSF fits data
        poscatalog : Catalog of positive detections.
        negcatalog : Catalog of negative detections.
             Sources extracted from the negative side
             of an image.
        snr_thresh : float, optional. Default is 100.
             Any source with 100 times the minimum SNR is
             considered a high SN source.
        local_thresh : float, optional. Default is 0.6.
             Sources with local variance greater than
             0.6 * negative noise are considered as 
             sources of high local variance.
        local_region : integer, optional. Default is 10.
             A region to compute the local variance in
             beam sizes.
        psfcorr : integer, optional. Default is 2.
             Data size to correlate. In beam sizes.
        high_corr_thresh :  float, optional. Default is 0.5.
             Correlation threshold. Sources of high correlation
             with the PSF have correlation > the specified.
        negdetec_region :  float, optional. Default is 10.
             Region to lookup for negative detections around
             a given source. In beam size.
        negative_thresh :  float, optional. Default is 6.
             Number of negative detections, N, threshold. Sources
             with number > N negatives around them are require direction
             dependent (DD) calibration solutions.
        phasecenter_excl_region :  float (in degrees), optional.
             A radius from the phase center (in beam sizes) to exclude 
             in making final DD source selection.
        prefix : str, optional. Sets a prefix to the output directory.
        loglevel :  int, optional. Default 0. Python logging.
        0, 1, 2, 3 for info, debug, error and critical respectively.
        """

        # image, psf image, positive and negative catalogues
        self.loglevel = loglevel
        self.prefix = prefix
        self.log = utils.logger(self.loglevel, prefix=self.prefix)
        self.imagename = imagename

        self.psfname = psfname
        if not self.psfname:
            self.log.info("dE tagging will be made without the PSF correlation"
                          " note that this might affect the results.")

        self.poscatalog = poscatalog
        self.negcatalog = negcatalog
        
        # reading the imagename data
        self.imagedata, self.wcs, self.header, self.pixsize =\
                          utils.reshape_data(self.imagename, prefix=self.prefix)
        self.log.info(" Loading image data")

        # computing the noise
        self.noise = utils.negative_noise(self.imagedata)
        self.log.info(" The negative noise of an image is %e Jy/beam"%
                       self.noise)

        # tags
        self.snr_tag = "snr"
        self.local_tag = "high_var"
        self.corr_tag = "high_corr"
        self.dd_tag = "dE"

        # thresholds
        self.snr_thresh = snr_thresh
        self.localthresh = local_thresh
        self.corrthresh = high_corr_thresh
        self.negthresh = negatives_thresh
       
        
        #regions
        self.psfcorr_region = psfcorr_region
        self.local_region = local_region
        self.phaserad = phasecenter_excl_radius # radius from the phase center
        self.negregion =  negdetec_region # region to look for negatives
        
        # central ra, dec, beam major axes
        self. ra0 = numpy.deg2rad(self.header["CRVAL1"])
        self.dec0 = numpy.deg2rad(self.header["CRVAL2"])
        self.bmaj = self.header['BMAJ'] # in degrees

        # Models
        self.pmodel = Tigger.load(self.poscatalog, verbose=self.loglevel)
        self.nmodel = Tigger.load(self.negcatalog, verbose=self.loglevel)
        self.r2d = 180.0/math.pi
        self.d2r = math.pi/180.0



    def number_negatives(self, source):
        

        #sources = filter(lambda src: src.getTag(tag), psources)

        tolerance = self.negregion * self.bmaj * self.d2r

        if self.phaserad:
            radius = numpy.deg2rad(self.phaserad * self.bmaj)
            
        
        ra, dec = source.pos.ra, source.pos.dec # in radians
        within = self.nmodel.getSourcesNear(ra, dec, tolerance)    
        if len(within) >= self.negthresh:
            if self.phaserad:
                if dist(self.ra0, dec, ra, self.dec0)[0] > radius: 
                        source.setTag(self.dd_tag, True)
            else:
                source.setTag(self.dd_tag, True)


    def source_selection(self):
        
        sources = self.pmodel.sources
            
        if self.noise == 0:
            self.log.error(" Division by 0. Aborting")

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
                            
        self.pmodel.save(self.poscatalog)
        return self.poscatalog
