# Makhuduga Lerato Sebokolodi <mll.sebokolodi#gmail.com>
# Computes the reliability of sources

import numpy
import pylab
from scipy import stats
import utils
import math
from astLib.astWCS import WCS
from collections import Counter
from scipy.interpolate import griddata
import pyfits


class KernelDensityEstimate(stats.gaussian_kde):
    # knicked from SoFiA
    def __init__(self, dataset, covariance):

        """ Makes a kernel density field of dataset with covariance matrix

        dataset: the data
        covariance: the covariance matrix
        """
        self.covariance = covariance
        stats.gaussian_kde.__init__(self, dataset)

    def _compute_covariance(self):
        self.inv_cov = pylab.linalg.inv(self.covariance)
        self._norm_factor = numpy.sqrt(
            pylab.linalg.det(2 * numpy.pi * self.covariance)) * self.n


class load(object):

    def __init__(self, positive_model, negative_model, image=None,
                 psf_image=None, do_local_variance=False,
                 do_psf_correlation=False, do_nearby_sources=False, local_var_step=10,
                 correlation_step=5, nearby_step=5, do_reliability_plots=True,
                 do_diag_cov=False, plot_prefix=None):

        """ Computes the reliability of the sources in positive model.

        positive_model: Takes in the positive model (either open or it will open it inside the script)
        negative_model: Takes in the model containing detections from the negative pixeled image
        image: If you have an image you can provide it. If not, no optional parameters will be added.
        psf_image: this is either the PSF  image



        """
        self.positive_model = positive_model
        self.negative_model = negative_model
        self.image = image
        self.psf_image = psf_image
        self.do_local = do_local_variance
        self.do_corr = do_psf_correlation
        self.do_near = do_nearby_sources
        self.near_step = nearby_step
        self.local_step = local_var_step
        self.psf_step = correlation_step
        self.do_plots = do_reliability_plots
        self.do_diag_cov = do_diag_cov
        self.prefix = plot_prefix
        utils.info(" >>>> Running reliability estimation script <<<<<")

        # if the image is not provided then no local variance or correlation
        # factor computation can be made.
        if self.image:
            with pyfits.open(self.image) as hdu:
                self.image_header = hdu[0].header
                self.image_data = utils.slice_image_data(hdu[0].data)
                self.beam_major = self.image_header["bmaj"]
                self.beam_minor = self.image_header["bmin"]
                self.wcs = WCS(self.image_header, mode="pyfits")
            if self.psf_image:
                with pyfits.open(self.psf_image) as psf:
                    self.psf_data = utils.slice_image_data(psf[0].data)
                    self.psf_header = psf[0].header
            else:
                utils.info("No PSF provided. Thus, no psf correlation will be made. "
                           "Provide image and psf image to add this as a parameter. ")
                self.do_corr = False
        else:
            self.do_local, self.do_corr, self.do_near = False, False, False
            self.image_data = None
            utils.info("No image data, thus, no extra parameter to be computed. ")

    def parameters(self, model, imagedata=None):

        labels = dict(size=(0, "Log$_{10}$( Source area [arcsec] )"),
                      peak=(1, "Log$_{10}$( Peak flux [Jy] )"),
                      total=(2, "Log$_{10}$( Total flux [Jy] )"))
        if self.do_corr:
            labels.update({"correlation_factor": (len(labels),
                          "Log$_{10}$ (Correlation Factor)")})
        if self.do_local:
            labels.update({"local_variance": (len(labels),
                           "Log$_{10}$(Local Variance)")})
        if self.do_near:
            labels.update({"nearby_sources": (len(labels),
                            "Log$_{10}$(Nearby Detections)")})

        nsrc = len(model.sources)
        out = numpy.zeros([nsrc, len(labels)])
        if self.do_corr:
            c0 = abs(self.psf_header["crpix2"])
            psf_region = self.psf_data[abs(c0-self.psf_step):c0+self.psf_step,
                         abs(c0-self.psf_step): c0 + self.psf_step].flatten()

        area, peak, total, local, corr, near = [], [], [], [], [], []
        for i, src in enumerate(model.sources):
            ra, dec = src.pos.ra, src.pos.dec # in radians
            if self.do_near:
                _near = model.getSourcesNear(ra, dec, self.near_step *
                                             (self.beam_major/3600.0))
                near.append(_near)
                src.setAttribute("nearby_sources", _near)

            ex = numpy.rad2deg(src.get_attr("_pybdsm_Maj"))
            ey = numpy.rad2deg(src.get_attr("_pybdsm_Min"))
            area_second = math.pi * ex * ey * pow(3600, 2) # converting from degrees to arc-second
            integrated_flux = src.flux.I
            peak_flux = src.getTag("_pybdsm_Peak_flux")
            ra_d, dec_d = map(lambda a: numpy.rad2deg(a), [ra, dec])
            if self.image:
                pix_ra, pix_dec = [self.wcs.wcs2pix(*(ra_d, dec_d))][0]
            area.append(area_second), peak.append(peak_flux), total.append(integrated_flux)

            if self.do_local:
                local_step = self.local_step * self.beam_major/\
                             float(abs(self.image_header['cdelt1']))
                local_data_region = self.image_data[abs(pix_dec-local_step):
                                    pix_dec+local_step, abs(pix_ra-local_step):
                                    pix_ra+local_step].flatten()
                # considering only the negative pixels.
                source_local = abs(local_data_region[local_data_region < 0]).std()
                local.append(source_local)
                src.setAttribute("local_variance", source_local)

            if self.do_corr:
                data_region = imagedata[abs(pix_dec - self.psf_step): pix_dec + self.psf_step,
                              abs(pix_ra - self.psf_step): pix_ra +self.psf_step].flatten()
                normalise_data = (data_region - data_region.min())/(data_region.max() -
                                                                    data_region.min())
                correlation = numpy.corrcoef((normalise_data, psf_region))
                corr_factor = ((numpy.diag((numpy.rot90(correlation))**2).sum())**0.5)/2**0.5
                corr.append(corr_factor)
                src.setAttribute("correlation_factor", corr_factor)

        for i, srs in enumerate(model.sources):
            if self.do_corr and self.do_local and self.do_near:
                out[i,...] = area[i], peak[i], total[i], corr[i], local[i], near[i]
            elif self.do_corr and self.do_local and not self.do_near:
                out[i,...] = area[i], peak[i], total[i], corr[i], local[i]
            elif self.do_corr and not self.do_local and self.do_near:
                out[i,...] = area[i], peak[i], total[i], corr[i], near[i]
            elif not self.do_corr and self.do_local and self.do_near:
                out[i,...] = area[i], peak[i], total[i], local[i], near[i]
            elif self.do_corr and not self.do_local and not self.do_near:
                out[i,...] = area[i], peak[i], total[i], corr[i]
            elif not self.do_corr and self.do_local and not self.do_near:
                out[i,...] = area[i], peak[i], total[i], local[i]
            elif not self.do_corr and not self.do_local and self.do_near:
                out[i,...] = area[i], peak[i], total[i], near[i]
            else:
                out[i,...] = area[i], peak[i], total[i]

        # removes the rows of 0s for the reliability
        remove_zeros = (out == 0).sum(1)
        output = out[remove_zeros <=0, :]
        return model, numpy.log10(output), labels

    def reliability_plots(self, positive, negative,
                          covariance, labels, png_out):
        utils.info(" >>>>>>>> Making the reliability plots <<<<<<<<<<")
        plots = []
        nplanes = len(labels)
        for i, label_i in enumerate(labels):
            for label_j in labels.keys()[i + 1:]:
                i = labels[label_i][0]
                j = labels[label_j][0]
                plots.append([i, j, label_i, label_j])
        npos, nneg = len(positive), len(negative)
        row  = 2
        column = 3
        pylab.figure(figsize=[20,10])
        cov = numpy.zeros([2, 2])
        for counter, (i, j, x, y) in enumerate(plots):
            pylab.subplot(int(row), int(column), counter + 1)
            a, b = negative[:, i], negative[:, j]
            c, d = positive[:, i], positive[:, j]
            cov[0, 0] = covariance[i, i]
            cov[0, 1] = covariance[i, j]
            cov[1, 0] = covariance[j, i]
            cov[1, 1] = covariance[j, j]

            ncov = KernelDensityEstimate(numpy.array([a, b]), cov)
            pn = ncov(numpy.array([a, b])) * nneg
            ac = numpy.concatenate((a, c))
            bd = numpy.concatenate((b, d))
            xi = numpy.linspace(ac.min(), ac.max(), 100)
            yi = numpy.linspace(bd.min(), bd.max(), 100)
            zzz = griddata((a,b), values=pn,
                           xi=(xi[None, :], yi[:, None]),
                           method='cubic')

            pylab.contour(xi, yi, zzz, 10, linewidths=2, colors='r')
            pylab.scatter(positive[:, i], positive[:, j], marker='o', c='k')
            pylab.tick_params(axis='x')
            pylab.tick_params(axis='y')
            pylab.xlabel(labels[x][1])
            pylab.ylabel(labels[y][1])
            pylab.xlim(ac.min(), ac.max())
            pylab.ylim(bd.min(), bd.max())
            pylab.grid()
        utils.info("Saving plot as %s"%png_out)
        pylab.savefig(png_out)

    def check_degeneracy(self, model_one, model_two, features):
        "Removes parameters columns which are degenerate"
        remove = []
        for x, key in enumerate(features.keys()):
            read_list = model_one[x, :]
            most_counter = Counter(read_list)
            mode, times = most_counter.most_common(1)[0]
            frac_rate = times/float(len(read_list))
            if frac_rate > 0.70:
                utils.info("Degenerate feature=%s is removed from reliablity estimation"%key)
                remove.append(x)
                del features[key]
            # delete the degenerate axis from it's model and the second model.
            model_one =  numpy.delete(model_one, remove, axis=0)
            model_two = numpy.delete(model_two, remove, axis=0)
        return features, model_one, model_two

    def get_reliability(self):
        # This removes sources at the edges.
        if self.do_local or self.do_corr:
            loc, cor = map(lambda a: (a * self.beam_major)/float(abs(self.image_header["cdelt1"])),
                           [self.local_step, self.psf_step])
            step = max([loc, cor])
            self.positive_model = utils.remove_sources_edges(
                self.positive_model, hdr=self.image_header, step=step)
            self.negative_model = utils.remove_sources_edges(
                   self.negative_model, hdr=self.image_header, step=step)

        pmodel, positive, labels = self.parameters(self.positive_model, self.image_data)
        negative_pixels = -1 * self.image_data if self.image_data is not None else None
        nmodel, negative, labels = self.parameters(self.negative_model, negative_pixels)
        utils.info("Checking degeneracies in the derived parameters. If parameter column"
                   " has 70 percent of the same value, the feature in question removed from the"
                   " reliability estimation.")
        features, positives, negatives = self.check_degeneracy(positive, negative, labels)

        bandwidth = []

        for plane in negatives.T:
            q3, q2 = numpy.percentile(plane, [75, 25])
            iqr = q3 - q2
            bandwidth.append(min([plane.std(), iqr/1.34]))
        nplanes = len(features)
        cov = numpy.zeros([nplanes, nplanes])
        nnsrc = len(negatives)
        npsrc = len(positives)
        utils.info("There are %d positive and %d negative detections to use for"
                       " reliability estimations."%(npsrc, nnsrc))
        # writing the diagonal part of the covariance matrix
        # using the normal reference rule (http://repec.org/esAUSM04/up.1603.1077410300.pdf)
        for i in range(nplanes):
            for j in range(nplanes):
                if i==j:
                    cov[i, j] = bandwidth[i] * ((4.0/((nplanes + 2) * nnsrc)) **
                                                   (1.0/(nplanes + 4.0)))
        if self.do_diag_cov:
            covariance = cov
            utils.info("Using a diagonal covariance matrix.")
        else:
            utils.info("Using a full covariance matrix.")
            diagonal = numpy.logical_not(numpy.diag(numpy.ones(nplanes))).astype(int)
            full_cov = numpy.cov(negatives.T)
            off_diagonal_cov = diagonal * full_cov
            covariance =  cov + off_diagonal_cov

        utils.info("The derived covariance matrix is: \n%s"%covariance)
        positive_density = KernelDensityEstimate(positives.T, covariance)
        negative_density = KernelDensityEstimate(negatives.T, covariance)
        nps = positive_density(positives.T) * npsrc
        nns = negative_density(positives.T) * nnsrc

        # computing the reliability
        reliability = (nps - nns)/ nps
        for src, rel in zip(self.positive_model.sources, reliability):
            src.setAttribute("rel", rel)
        utils.info("Saving the reliability values....")

        if self.do_plots:
            png_ = self.prefix + "_rel_plots.png" if self.prefix \
                else "reliability_plots.png"
            self.reliability_plots(positives, negatives, covariance,
                               labels=features, png_out= png_)

        return self.positive_model, self.negative_model
    # save the tags outside




























