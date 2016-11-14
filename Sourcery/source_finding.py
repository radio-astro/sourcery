from lofar import bdsm
import utils


def find_sources(image, output, catalog_type="gaul",
                 format="ascii", **kw):
    
    """This calls pybdsm to find sources on an image

    image: an input image to run source finder.

    kw: excepts all PyBDSM extra parameters. For more
    information see PyBDSM manual.

    """

    if not isinstance(image, str):
        utils.abort("Image provided is not a string. %s"%image)

    ext = utils.file_ext(image)
    img = bdsm.process_image(image, **kw)

    img.write_catalog(outfile=output, catalog_type=catalog_type,
    format=format, clobber=True)

    return output











