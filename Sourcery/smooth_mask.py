import pyfits
import numpy
from scipy.ndimage import filters
import math
import utils


def smooth_and_mask(image, output_mask, smooth=True,
                    scales=None, thresh=2, noise=None):

    """ This function smoothness the image data at different scales which a user can provide.
        Defaults are: [0.1, 0.4, 0.5, 1.0, 1.5, 2, 2.5]

        image: the image to smooth and/or mask
        output_mask: the mask data
        smooth: if true, smoothing is applied before masking
        scales: are the sigma scales to use during Gaussian filtering
        thresh: a factor above the noise to mask, default is 2.
        noise: pass noise estimate if you know it else it will be computed inside program.
    """

    with pyfits.open(image) as hdu:
        data = hdu[0].data
        hdr = hdu[0].header

    ndim = hdr["naxis"]
    kernel = numpy.ones(ndim, dtype=int).tolist()
    mask = numpy.ones(data.shape)
    masked = 0

    noise = noise or utils.average_negative_noise(data)[0]
    threshold = noise * thresh
    if smooth:
        scale = scales or [0.1, 0.3, 0.5, 1.0, 1.5, 2, 2.5]
        utils.info("Smoothing using scales = %s"%scale)
        emin, emaj = hdr["bmin"], hdr["bmaj"]
        cell = abs(hdr["cdelt1"])
        beam = math.sqrt( emin * emaj )/float(cell)
        smooth = None
        for sc in scale:
            kk = sc * beam
            kernel[-1] = kk
            kernel[-2] = kk
            smooth_data = filters.gaussian_filter(data if smooth is None else smooth, kernel)
            mask *= smooth_data > threshold
            masked += mask # adding smoothed masked data

        # converting all the masked data to 0s and 1s
        masked[masked > 1] = 1
        masked_data = masked * data
    else:
        utils.info("Creating a mask without smoothing.")
        mask = data > threshold
        masked_data = mask * data

    hdu[0].data = masked_data
    hdu.writeto(output_mask, clobber=True) # clobber true -- replaces any existing mask with similar name
    return output_mask

