# Determing sources that require direction dependent calibration solutions
# Lerato Sebokolodi <mll.sebokolodi@gmail.com>


import matplotlib
matplotlib.use('Agg')
import Tigger
from Tigger.Coordinates import angular_dist_pos_angle as dist
import utilss
import numpy


class Parent(object):


    def __init__(self, imagename, psfname, poscatalog, negcatalog,
                 snr_thresh=0.5, snr_tag='snr', local_thresh=0.9,
                 local_region=10, high_local_tag='high_variance',
                 psfcorr_region=2, high_corr_thresh=0.4,  high_corr_tag='cf',
                 negradius=10, negatives_thresh=6, excl_sources_radius=None,
                 dd_tag='dE', loglevel=0):


        """ Determines sources that require direction-dependent (DD)
        calibration solutions.

        imagename: Fits data
        psfname : PSF fits data
        poscatalog : Catalog of positive detections.
        negcatalog : Catalog of negative detections.
             Sources extracted from the negative side
             of an image.
        snr_thresh : float, optional. Default is 0.5.
              Any sources 0.5 times the minimum SNR is
              considered a high SN source.
        snr_tag : str, optional. Default is 'snr'.
             If a tag specified all sources of high SNR
             will tagged as such.

        local_thresh : float, optional. Default is 0.9.
             Sources with local variance greater than
             0.9* negative noise are considered as 
             high local variance sources.
        local_region : integer, optional. Default is 10.
             A region to compute the local variance, e.g.,
             local_size = 10, means that 10 * beam size 
             region will be used to compute the local
             variance.
        high_local_tag : string, optional. Default is 'high_variance'.
             A tag given to sources of high variance.

        psfcorr : integer, optional. Default is 2.
             Specifies the size data to correlated, e.g
             2 implies 2 times the beam size data to
             correlate.
        high_corr_thresh :  float, optional. Default is 0.4.
             Sources with correlation with the instrument's
             PSF larger than 0.4 are considered as high correlation.
        high_corr_tag :  string, optional. Default is 'cf'.
             Sources of high correlation are tagged as such.
        negradius :  float, optional. Default is 10 in degrees.
             Sets a tolerance radius to look up for negative
             detections around a given source, e.g.,
             for negradius = 10 then search for negatives
             around a single position detection at a 
             radius = 10 * beam major axes i.e 10 times beam
             sizes.
        negative_thresh :  float, optional. Default is 6.
             If a source has N number of negative detections 
             around it such that N > 6 then that source requires
             the direction-dependent calibration.
        excl_sources_radius :  float (in degrees), optional. Default is None.
             Defines a radius to sources to exclude when making
             DD source selections, e.g., given a value 5 then sources
             with distance < 5 * beam size from the image center are
             excluded when making the final dd_tag (see below) tagging.
        dd_tag : str, optional. Default is 'dE'.
             Sources that requires DD calibration are given
             this tag. 
                
        """

        self.imagename = imagename
        self.psfname = psfname
        self.poscatalog = poscatalog
        self.negcatalog = negcatalog
        self.loglevel = loglevel
        self.log = utilss.logger(self.loglevel)


        self.imagedata, self.wcs, self.header, self.pixsize=\
             utilss.reshape_data(self.imagename)
        self.log.info("Loading image data")

        self.noise = utilss.negative_noise(self.imagedata)
        self.log.info("The negative noise of an image is %e"%
                     self.noise)

        self.snr_tag = snr_tag
        self.snr_thresh = snr_thresh

        self.local_thresh = local_thresh
        self.high_local_tag = high_local_tag
        self.local_region = local_region

        self.high_corr_thresh = high_corr_thresh
        self.high_corr_tag = high_corr_tag
        self.psfcorr_region = psfcorr_region

        self.negatives_thresh = negatives_thresh
        self.excl_sources_radius = excl_sources_radius
        self.dd_tag = dd_tag
        

        self. ra0 = numpy.deg2rad(self.header["CRVAL1"])
        self.dec0 = numpy.deg2rad(self.header["CRVAL2"])
        self.bmaj_rad = numpy.deg2rad(self.header['BMAJ'])
        self.bmaj_pix = self.header['BMAJ']/ self.pixsize



        def signal_to_noise(self):

            model = Tigger.load(self.poscatalog)
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

            pmodel = Tigger.load(self.poscatalog)
            nmodel = Tigger.load(self.negcatalog)

            psources = pmodel.sources

            sources = filter(lambda src: src.getTag(self.high_corr_tag),
                            psources)
      
            tolerance = numpy.deg2rad(self.negradius) * self.bmaj_rad
            
            if self.excl_sources_radius:
                radius = numpy.deg2rad(self.excl_sources_radius) * self.bmaj_rad 
            for srs in sources:
                ra, dec = srs.pos.ra, srs.pos.dec
                within = nmodel.getSourcesNear(ra, dec, tolerance)
                if len(within) > self.negatives_thresh:
                    if self.excl_sources_radius:
                        if dist(ra, self.ra0, dec, self.dec0)[0] > radius : 
                            srs.setTag(self.dd_tag, True)
                    else:
                            srs.setTag(self.dd_tag, True)
            pmodel.save(self.poscatalog)

        
        def source_selection(self):
             
            signal_to_noise()

            utilss.local_variance(self.imagedata, self.header,
                  self.poscatalog, self.wcs, self.pixsize, tag=self.snr_tag,
                  local_region=self.local_region, noise=self.noise,
                  highvariance_factor=self.local_thresh, high_local_tag=
                  self.high_local_tag, neg_side=True, do_alltag=False,
                  do_high_loc=True)

            utilss.psf_image_correlation(catalog=self.poscatalog,
                  psfimage=self.psfname, imagedata=self.imagedata,
                  header=self.header, wcs=self.wcs, pixelsize=
                  self.pixsize, corr_region=self.psfcorr_region,
                  thresh=self.high_corr_thresh, tags=self.high_local_tag,
                  coefftag=self.high_corr_tag, do_alltag=False, do_high=True)
             
            number_negatives()


            


