import time
import numpy
import subprocess
import os
import sys
from astLib.astWCS import WCS
import tempfile
import pyfits
import datetime

def info(string):
    t = "%d/%d/%d %dh:%dm:%ds"%(time.localtime()[:6])
    print "%s ##INFO: %s"%(t, string)


def warn(string):
    t = "%d/%d/%d %dh:%dm:%ds"%(time.localtime()[:6])
    print "%s ##WARNING: %s"%(t, string)


def abort(string, exception=None):
    t = "%d/%d/%d %dh:%dm:%ds"%(time.localtime()[:6])
    exception = exception or SystemExit
    raise exception("%s ##ABORTING: %s"%(t, string))


def file_ext(filename):
    ext = filename.split(".")[-1]
    return ext


def slice_image_data(data):

    """ returns first image slice of data """
    imslice = numpy.zeros(data.ndim, dtype=int).tolist()
    imslice[-1] = slice(None)
    imslice[-2] = slice(None)

    return data[imslice]


def average_negative_noise(imagedata):

    """Computes the noise of an image from the negative pixels only"""
    if isinstance(imagedata, str):
        with pyfits.open(imagedata) as hdu:
            data = hdu[0].data
    else:
        data = imagedata
    data2d = slice_image_data(data)
    negative = data2d[data2d < 0].flatten()
    noise_std = numpy.concatenate([negative, -negative]).std()
    noise_mean = numpy.mean(abs(negative.flatten() - numpy.mean(negative.flatten())))
    info("Noise compute successful: Negative noise std = %e  and noise mean deviation = %e"
         %(noise_std, noise_mean))

    return noise_std, noise_mean


def remove_sources_edges(model, hdr, step):
    """Removes sources at the edges from the model.

    model: is model of sources
    hdr: header of the image
    step: the step size in pixels.

    """
    naxis = hdr["naxis1"]
    xmax = naxis - step
    xmin = step
    wcs = WCS(hdr, mode="pyfits")
    for i, src in enumerate(model.sources):
        ra, dec = map(lambda a: numpy.rad2deg(a), [src.pos.ra, src.pos.dec])
        rap, decp = [wcs.wcs2pix(*(ra, dec))][0]
        if rap < xmin or decp < xmin or rap > xmax or decp > xmax:
            model.sources.remove(src)
    return model


def get_prefix(prefix, imagename, outdir):

    if prefix is None:
        # base it on the name
        prefix = os.path.basename(imagename.split(",")[0]).split(".")[:-1]
        prefix = ".".join(prefix)
    # set the out directory
    outdir = outdir or prefix + "_" + os.path.join(datetime.datetime.now(). \
                                                   strftime("%Y-%m-%d-%H"))
    outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    prefix = outdir + "/" + prefix
    return prefix



def make_inverted_image(image, prefix, ext='fits', save=False):
    info(" Making the inverted image. ")
    if save:
        negative_image = prefix + "_negative." + ext
    else:
        negative_image = tempfile.NamedTemporaryFile(suffix="." + ext, dir=".")
        negative_image.flush()
        negative_image = negative_image.name
    with pyfits.open(image) as hdu:
        image_data = slice_image_data(hdu[0].data)
    pyfits.writeto(negative_image, -image_data, hdu[0].header, clobber=True)
    return negative_image


def xrun(command, options, log=None):
    """
        Run something on command line.
        Example: _run("ls", ["-lrt", "../"])
    """

    options = map(str, options)

    cmd = " ".join([command] + options)

    if log:
        log.info("Running: %s" % cmd)
    else:
        print('running: %s' % cmd)

    process = subprocess.Popen(cmd,
                               stderr=subprocess.PIPE if not isinstance(sys.stderr, file) else sys.stderr,
                               stdout=subprocess.PIPE if not isinstance(sys.stdout, file) else sys.stdout,
                               shell=True)

    if process.stdout or process.stderr:

        out, err = process.comunicate()
        sys.stdout.write(out)
        sys.stderr.write(err)
        return out
    else:
        process.wait()
    if process.returncode:
        raise SystemError('%s: returns errr code %d' % (command, process.returncode))























