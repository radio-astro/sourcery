# Lerato Sebokolodi <mll.sebokolodi@gmail.com>

import sys
from argparse import ArgumentParser
import source_finding
import utils
import compute_reliability
import os
import simplejson as json
import smooth_mask
import Tigger
import tempfile


# set sourcery directory
_path = os.path.realpath(__file__)
_path = os.path.dirname(_path)
execfile("%s/__init__.py"%_path)


def main():

    for i, arg in enumerate(sys.argv):
        if (arg[0] == '') and arg[1].isdigit(): sys.argv[i] = '' + arg

    parser = ArgumentParser(description="Computes Reliablities and Tags Sources that require"
                           " Direction-Dependent correction")
    add = parser.add_argument

    add("-v", "--version", action="version", version=
        "{:s} version {:s}".format(parser.prog, __version__))

    add("-i", "--image", dest="image", help="Fits image to extract sources from.",
        type=str, default=None)

    add("-p", "--psf-image", dest="psf_image", default=None, type=str)

    add("-pc", "--pos-catalog", dest="positive_catalog", default=None, type=str,
        help=" Takes in a model of detections obtained from the positive image. "
             "Current format supported: Tigger'.lsm.html'")

    add("-nc", "--neg-catalog", dest="negative_catalog", default=None, type=str,
        help=" Model containing detections from the negative pixels (in Tigger format)")

    add("-pref", "-dir-prefix", dest="prefix", default=None,
        help=" Prefix for output products. System default will be generated.")

    add("-od", "--output-dir", dest="outdir", default=None,
        help=" Output products will be dumped here. System Default will be generated.")

    add("-adl", "--do-local", dest="add_localvar", action="store_true",
        default=False, help="Include local variance for reliability estimations."
                            "Default is False")

    add("-apsf", "--do-psf", dest="add_psf_corr", action="store_true",
        default=False, help="Add PSF correlation for reliability estimations."
                            "Default is False")

    add("-ans", "--do-nearby-sources", dest="add_near", action="store_true",
        default=False, help="Add number of nearby sources as a feature. Default is False.")

    add("-lr", "--local-region", dest="local_region", default=10,
        help="The size of image pixels to compute the local variance. E.g"
             " if 5 is given, then the region of 5 * PSF sizes  centred about"
             " a source will be used.  Default is 10.")

    add("-pr", "--psf-region", dest="psf_region", default=5,
        help= "The size of the image pixels to correlate with the PSF, see -lr."
              " Default is 5.")

    add("-nr", "--nearby-region", dest="nearby_region", default=10,
        help="The region to look up for source, see -lr. Default is 10.")

    add("-drp", "--do-rel-plots", dest="do_rel_plots", default=True,
        help="Make reliability plots. Default is True.", action="store_false")

    add("-dg", "--diag-cov", dest="use_diag_cov", default=False, action="store_true",
        help="If specified, a diagonal covariance matrix will be used to estimate the density"
             " functions.")

    add("--to-pybdsm", dest="to_pybdsm", action="append",
           help=" PyBDSM process_image option, e.g. thresh_pix=1,thresh_pix=2")

    add("-cat", "--catalog-type", dest="catalog_type", default="gaul",
           help="A catalog type (pybdsm module), e.g gaul, srl"
                " which corresponds to Gaussian and Source components. This option"
                " is available if source finding is performed within the program. Default"
                " is gaul.")

    add("-fm", "--catalog-format" ,dest="catalog_format", default="ascii",
           help="Format to save your catalog. The available options are "
                "ascii (default) and fits.")

    add("-sn", "--save-neg-image", dest="save_negative_image", default=True,
        help="Save the inverted image. Default is True.", action="store_false")

    add("-smt", "--smooth-mask", dest="smooth_mask", default=False,
        help=" Smooth the image using Gaussian kernels at different scales."
             " Then mask the pixels a factor c below the noise."
             " Default is False.", action="store_true")

    add("-mf", "--mask-thresh", dest="mask_thresh", default=2.0,
        help="This will mask all pixels in the smoothed image below "
             "this factor above the noise. Default is 2.")

    add('-jc', '--json-config', dest='config', default=None,
        help='Json config file : No default')

    args = parser.parse_args()
    if args.smooth_mask:
        smooth_fits = tempfile.NamedTemporaryFile(suffix="." + "fits", dir=".")
        smooth_fits.flush()
    if args.config:
        with open(args.config) as conf:
            jdict = json.load(conf)
            for key, val in jdict.items():
                if isinstance(val, unicode):
                    jdict[key] = str(val)
        # reading in the different dictionaries in the config file
        reliability_dict = jdict["reliability"]
        dde_tagging_dict = jdict["dd_tagging"]
        source_finding_dict = jdict["source_finder_opts"]

        outdir = jdict["outdir"] or args.outdir
        cat_type = jdict["catalog_type"] or args.catalog_type
        cat_format = jdict["catalog_format"] or args.catalog_format
        prefix = jdict["prefix"] or args.prefix
        positive_catalog = jdict["positive_catalog"] or args.positive_catalog
        negative_catalog = jdict["negative_catalog"] or args.negative_catalog
        images  = jdict["image"] or args.image
        psf_images = jdict["psf_image"] or args.psf_image

        if positive_catalog and negative_catalog:
            image = images.split(",")[0]
            psf_image = psf_images.split(",")[0] if psf_images else [None]
            prefix = utils.get_prefix(prefix, positive_catalog, outdir)
            pmodel = Tigger.load(positive_catalog)
            nmodel = Tigger.load(negative_catalog)
            reliability_dict["plot_prefix"] = prefix
            load_reliability = compute_reliability.load(
                positive_model=pmodel, negative_model=nmodel, image=image,
                psf_image=psf_image, **reliability_dict)
            positive_model, negative_model = load_reliability.get_reliability()
            positive_model.save(prefix + ".lsm.html")
            negative_model.save(prefix + "_negative.lsm.html")
            utils.info("Sourcery Completed successfully!!!")

        elif images:
            images = images.split(",")
            psfs = psf_images.split() if psf_images else [None]
            prefix = utils.get_prefix(prefix, images[0], outdir)
            if len(images) != len(psfs):
                psfs = [psfs[0]]*len(images)
            if args.smooth_mask or jdict["smooth_mask"]:
                smooth_fits = tempfile.NamedTemporaryFile(suffix="." + "fits", dir=".")
                smooth_fits.flush()
            for i, (image, psf) in enumerate(zip(images, psfs)):
                if len(images) > 1:
                    prefix = prefix + "-%04d"%i
                reliability_dict["plot_prefix"]  = prefix
                ext = ext = utils.file_ext(images[0])
                output = prefix + "." + cat_type

                if args.smooth_mask or jdict["smooth_mask"]:
                    noise = utils.average_negative_noise(image)[0]
                    smooth_mask.smooth_and_mask(image, smooth_fits.name,
                                                thresh=args.mask_thresh,
                                                noise=noise)
                    source_finding_dict["blank_limit"] = 1e-9
                    source_finding_dict["detection_image"] = smooth_fits.name
                positive_catalog = source_finding.find_sources(
                    image, output=output, catalog_type=cat_type,
                format=cat_format, **source_finding_dict)
                negative_image = utils.make_inverted_image(
                    image, prefix, ext=ext, save=args.save_negative_image)
                source_finding_dict['thresh_isl'] = int(source_finding_dict['thresh_isl']) - 1 if \
                    source_finding_dict.get('thresh_isl') else 2
                source_finding_dict['thresh_pix'] = int(source_finding_dict['thresh_pix']) - 1 \
                    if source_finding_dict.get('thresh_pix') else 4
                if args.smooth_mask or jdict["smooth_mask"]:
                    smooth_mask.smooth_and_mask(negative_image, smooth_fits.name,
                                                thresh=args.mask_thresh, noise=noise)
                    source_finding_dict["detection_image"] = smooth_fits.name
                negative_catalog = source_finding.find_sources(
                    negative_image, output.replace(
                        "."+cat_type, "_negative."+cat_type), **source_finding_dict)
                pmodel = Tigger.load(positive_catalog)
                nmodel = Tigger.load(negative_catalog)
                load_reliability = compute_reliability.load(
                    positive_model=pmodel, negative_model=nmodel,
                    image=image, psf_image=psf, **reliability_dict)
                positive_model, negative_model = load_reliability.get_reliability()
                positive_model.save(prefix + ".lsm.html")
                negative_model.save(prefix + "_negative.lsm.html")
            if args.smooth_mask or jdict["smooth_mask"]:
                smooth_fits.close()
        else:
            utils.abort("No images or models provided on the terminal or config file. Exiting.")

    else:
        # taking the pybdsm options from the command line.
        to_pybdsm = args.to_pybdsm[0].split(",") if args.to_pybdsm else None
        pybdsm_opts = dict([items.split("=") for items in to_pybdsm]) \
            if to_pybdsm else {}
        cat_type = args.catalog_type  # gaul, srl
        cat_format = args.catalog_format  # ascii, fits
        # checking the entries of pybdsm_opts
        for key, val in pybdsm_opts.iteritems():
            if not callable(val):
                try:
                    pybdsm_opts[key] = eval(key)
                except NameError:
                    pybdsm_opts[key] = val

        if args.positive_catalog and args.negative_catalog:
            utils.info("Using models %s and %s for the reliablity estimation " %
                       (args.positive_catalog, args.negative_catalog))
            prefix = utils.get_prefix(args.prefix, args.positive_catalog, args.outdir)
            pmodel = Tigger.load(args.positive_catalog)
            nmodel = Tigger.load(args.negative_catalog)
            load_reliability = compute_reliability.load(
                positive_model=pmodel, negative_model=nmodel, image=args.image,
                psf_image=args.psf_image, do_local_variance=args.add_localvar,
                do_psf_correlation=args.add_psf_corr, do_nearby_sources=args.add_near,
                local_var_step=args.local_region, correlation_step=args.psf_region,
                nearby_step=args.nearby_region, do_reliability_plots=args.do_rel_plots,
                do_diag_cov=args.use_diag_cov, plot_prefix=prefix)
            positive_model, negative_model = load_reliability.get_reliability()
            positive_model.save(prefix + ".lsm.html")
            negative_model.save(prefix + "_negative.lsm.html")
            utils.info("Sourcery Completed successfully!!!")
        elif args.image:
            images = args.image.split(",")
            psfs = args.psf_image.split(",") if args.psf_image else [None]
            ext = utils.file_ext(images[0])
            if len(images) != len(psfs):
                psfs = [psfs[0]]*len(images)

            utils.info("Running Source Finding on the image(s) provided. ")
            for i, (image, psf) in enumerate(zip(images, psfs)):
                prefix = utils.get_prefix(args.prefix, images[0], args.outdir)
                if len(images) > 1:
                    prefix = prefix + "-%04d"%i
                output = prefix + "." + cat_type
                if args.smooth_mask:
                    noise = utils.average_negative_noise(image)[0]
                    smooth_mask.smooth_and_mask(image, smooth_fits.name,
                                                thresh=args.mask_thresh, noise=noise)
                    pybdsm_opts["detection_image"] = smooth_fits.name
                    pybdsm_opts["blank_limit"] = 1e-9
                positive_catalog = source_finding.find_sources(
                    image, output=output, catalog_type=cat_type,
                    format=cat_format, **pybdsm_opts)
                negative_image = utils.make_inverted_image(image, prefix, ext=ext,
                                                 save=args.save_negative_image)
                pybdsm_opts['thresh_isl'] = int(pybdsm_opts['thresh_isl']) - 1 if\
                    pybdsm_opts.get('thresh_isl') else 2
                pybdsm_opts['thresh_pix'] = int(pybdsm_opts['thresh_pix']) - 1 \
                    if pybdsm_opts.get('thresh_pix') else 4
                if args.smooth_mask:
                    smooth_mask.smooth_and_mask(negative_image, smooth_fits.name,
                                                thresh=args.smooth_mask, noise=noise)
                    pybdsm_opts["detection_image"] = smooth_fits.name
                negative_catalog = source_finding.find_sources(
                    negative_image, output=output.replace("." + cat_type, "_negative."+cat_type),
                    catalog_type=cat_type, format=cat_format, **pybdsm_opts)
                pmodel = Tigger.load(positive_catalog)
                nmodel = Tigger.load(negative_catalog)
                load_reliability = compute_reliability.load(
                            positive_model=pmodel, negative_model=nmodel, image=image,
                            psf_image=psf, do_local_variance=args.add_localvar,
                            do_psf_correlation=args.add_psf_corr, do_nearby_sources=args.add_near,
                            local_var_step=args.local_region, correlation_step=args.psf_region,
                            nearby_step=args.nearby_region, do_reliability_plots=args.do_rel_plots,
                            do_diag_cov=args.use_diag_cov, plot_prefix=prefix)
                positive_model, negative_model = load_reliability.get_reliability()
                positive_model.save(prefix + ".lsm.html")
                negative_model.save(prefix + "_negative.lsm.html")
                utils.info("Sourcery Completed successfully!!!")
            if args.smooth_mask:
                smooth_fits.close()
        else:
            utils.abort("No image or models provided. Aborting. See help command.")









