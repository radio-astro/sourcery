# Determing sources that require direction dependent calibration solutions
# Lerato Sebokolodi <mll.sebokolodi@gmail.com>


import matplotlib
matplotlib.use('Agg')
import Tigger
from Tigger.Coordinates import angular_dist_pos_angle as dist
import utils
import numpy


class load(object):


    def __init__(self, imagename, psfname, poscatalog, negcatalog,
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
        self.imagename = imagename
        self.psfname = psfname
        self.poscatalog = poscatalog
        self.negcatalog = negcatalog
        self.loglevel = loglevel
        self.prefix = prefix
        self.log = utils.logger(self.loglevel, prefix=self.prefix)
        
        # reading the imagename data
        self.imagedata, self.wcs, self.header, self.pixsize =\
                          utils.reshape_data(self.imagename, prefix=self.prefix)
        self.log.info("Loading image data")

        # computing the noise
        self.noise = utils.negative_noise(self.imagedata)
        self.log.info("The negative noise of an image is %e"%
                       self.noise)

        # tags
        self.snr_tag = "snr"
        self.high_local_tag = "high_var"
        self.high_corr_tag = "high_corr"
        self.dd_tag = "dE"

        # thresholds
        self.snr_thresh = snr_thresh
        self.local_thresh = local_thresh
        self.high_corr_thresh = high_corr_thresh
        self.negatives_thresh = negatives_thresh
        
        #regions
        self.psfcorr_region = psfcorr_region
        self.local_region = local_region
        self.phasecenter_excl_radius = phasecenter_excl_radius
        self.negdetec_region =  negdetec_region
        
        # central ra, dec, beam major axes
        self. ra0 = numpy.deg2rad(self.header["CRVAL1"])
        self.dec0 = numpy.deg2rad(self.header["CRVAL2"])
        self.bmaj_deg = self.header['BMAJ'] # in degrees


    def signal_to_noise(self):
        
        model = Tigger.load(self.poscatalog, verbose=self.loglevel)
        sources = model.sources
            
        if self.noise == 0:
            self.log.error("Division by 0. Aborting")

        snr = [src.flux.I/self.noise for src in sources]
            
        thresh = self.snr_thresh * min(snr)
        n = 0
        for srs, s in zip(sources,snr):
            if s > thresh:
                srs.setTag(self.snr_tag, True)
                n += 1
        self.log.info("There are %d with high SNR"%n)
        model.save(self.poscatalog)


    def number_negatives(self):
        
        pmodel = Tigger.load(self.poscatalog, verbose=self.loglevel)
        nmodel = Tigger.load(self.negcatalog, verbose=self.loglevel)
        psources = pmodel.sources 
        sources = filter(lambda src: src.getTag(self.high_corr_tag), psources)

        tolerance = numpy.deg2rad(self.negdetec_region * self.bmaj_deg)

        if self.phasecenter_excl_radius:
            radius = numpy.deg2rad(self.phasecenter_excl_radius * self.bmaj_deg)
            
        for srs in sources:
            ra, dec = srs.pos.ra, srs.pos.dec
            within = nmodel.getSourcesNear(ra, dec, tolerance)    
            if len(within) >= self.negatives_thresh:
                if self.phasecenter_excl_radius:
                    if dist( self.ra0, dec, ra, self.dec0)[0] > radius: 
                        srs.setTag(self.dd_tag, True)
                else:
                    srs.setTag(self.dd_tag, True)
        pmodel.save(self.poscatalog)

        
    def source_selection(self):
             
        # signal-to-noise ratio
        self.signal_to_noise()
        # local variance
        utils.local_variance(
             self.imagedata, self.header, self.poscatalog, self.wcs,
             self.pixsize, tag=self.snr_tag,local_region=self.local_region,
             noise=self.noise, highvariance_factor= self.local_thresh,
             high_local_tag=self.high_local_tag, neg_side=True,
             setatr=False, prefix=self.prefix, do_high_loc=True)
        # correlation
        utils.psf_image_correlation(
            catalog=self.poscatalog, psfimage=self.psfname,
            imagedata=self.imagedata, header=self.header,
            wcs=self.wcs, pixelsize=self.pixsize, corr_region=
            self.psfcorr_region, thresh=self.high_corr_thresh,
            tags=self.high_local_tag, coefftag=self.high_corr_tag,
            setatr=False, do_high=True, prefix=self.prefix)
        # number of negative detections
        self.number_negatives()

        return self.poscatalog, self.negcatalog
