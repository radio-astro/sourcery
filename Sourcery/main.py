# Reliability estimation and direction-dependent source selection
## Lerato Sebokolodi <mll.sebokolodi@gmail.com>

import sys
from argparse import ArgumentParser
import reliabilityestimates as rel
import direcdepen as dd
import os
import datetime
import json




# set simms directory
_path = os.path.realpath(__file__)
_path = os.path.dirname(_path)
execfile("%s/__init__.py"%_path)


def main():

    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

    parser = ArgumentParser(description="Reliability estimator and"
                           " direction-dependent sources selection tool."
                           " M. L. L Sebokolodi <mll.sebokolodi@gmail.com>")

    add = parser.add_argument
    add("-v","--version", action="version",version=
        "{:s} version {:s}".format(parser.prog, __version__))
   
    add("-i", "--image", dest="image",
        help="FITS image name")

    add("-p", "--psfimage", dest="psf", default=None,
        help="Point Spread Function (PSF) Fits image name. Default=None")

    add("-s", "--source-finder", dest="source_finder", default="pybdsm",
        help="Source finder name. Current available"
        " sourcefinders: pybdsm. Default='pybdsm'")

    add("--to-pybdsm", dest="to_pybdsm", action="append",
           help="PyBDSM process_image options" 
           " [may be specified multiple times]. E.g thresh_pix=1 etc.")

    add("-log", "--log-level",dest="log_level", type=int, default=0,
        help="Python logging module. This ranges between"
        " 0-3, where 0, 1, 3, 4 are for info, debug, error and " 
        " critical respectively. Default is 0")

    add("-od", "--output-dir", dest="outdir", default=None, 
        help="Output products will be dumped here. System Default will be generated")

    add("-pref", "-dir-prefix", dest="prefix", default=None, 
        help="Prefix for output products. System default will be generated")

    add("-apsf", "--add-psfcorr", dest="add_psfcorr", action="store_true",
        default=False, help="Do and add correlation of the sources with the"
        "PSF as an extra source parameter for reliability estimations.  The psf"
        " name must be provided. Default is False. To set true add -apsf")

    add("-alv", "--add-localvar", dest="add_locvar", action="store_true",
        default=False, help=" Include local variance as an extra source"
        "parameter. See -apsf. Similar to -apsf.")

    add("-dmp", "--do-relplots", dest="do_relplots", action="store_false",
        default=True, help="Make reliability density plot. Default is True."
        " To disable add -dmp on the command line.")
 
    add("-pcr", "--psfcorr-region",  dest="psfcorr_region", type=int,
        default=5, help="Data size to correlate, given in beam sizes."
        " Default is 5 beam sizes.")
  
    add("-lr", "--local-region", dest="local_region", default=10,
        help="Data size to load the local variance, in beam sizes."
        " Default is 10 beam sizes.")
     
    add("-rel_rm", "--rel-sources-excl", dest="rel_src_excl", 
        action="append", default=None, help="Remove sources in a given"
        "region to exclude in reliability estimations."
        "E.g ra, dec, radius (in degrees). For more than"
        " one region: ra1,dec1,radius1:ra2,dec2,radius2. Default is None.")

    add("-ps", "--positive-smooth", dest="pos_smooth", type=float, default=1.6, 
        help="Data smoothing to eliminate noise/data peaks given by -ps * noise."
        " This is for positive side of an Fits image. Default is 1.6.")

    add("-ns", "--negative-smooth", dest="neg_smooth", type=float, default=0.8, 
        help="This is similar to -ps above but for negative side"
        "of an image. Default is 0.8.")

    add('-pisl', "--thresh-isl", dest="thresh_isl", type=float, default=3,
        help="Threshold for the island boundary in number of sigma above"
        " the mean. Determines extent of island used for fitting [pybdsm]."
        " For positive pixels. Default is 3.")

    add("-ppix", "--thresh-pix", dest="thresh_pix", type=float, default=5,
        help="Source detection threshold: threshold for the island peak"
        " in number of sigma above the mean. For positive pixels."
        " Default is 5.")
 
    add("-nisl", "--negthreshold-isl", dest="neg_thresh_isl",
        type=float, default=3, help="Similar to -pisl but applied"
        "to the negative pixels. Default is 3.")

    add("-npix", "--negthreshold-pix", dest="neg_thresh_pix",
        type=float, default=5, help="Similar to -ppix but applied"
        "for negative pixels. Default is 5.")

    add("-snr_thr", "--snr-threshold", dest="snr_thresh", type=float,
        default=80, help="Signal-to-noise threshold with reference"
        "to the minimum source flux in an image, e.g 80* min(snr), sources with"
        " snr > than this are referred to high SNR sources."
        " Default is 80.")

    add("-loc_thr", "--localvar-threshold", dest="locvar_thresh",
        type=float, default=0.9, help="Local variance threshold."
        "For -loc-thr of 0.9 means that"
        " sources with local variance > 0.9 * negative noise"
        " are considered sources of high local variance."
        " Default is 0.9")

    add("-pc_thr", "--psfcorr-threshold", dest="psfcorr_thresh", 
        type=float, default=0.5, help="Correlation factor threshold."
        "Sources with threshold larger than the specified are"
        "considered sources of high correlation. Default is 0.5.")    

    add("-nneg", "--num-negatives",dest="num_negatives", type=float, default=8,
        help="Number of negative detections around a given source."
        " If N > threshold then sources will be considered as requiring"
        " direction-dependent calibratio solutions. Default is 8.")

    add("-nrgn", "--neg-region", dest="neg_region", type=float, default=10,
        help="The size of a region around a sources to lookup for"
        " negatives. In Beams in sizes. Default is 10.") 

    add("-nphrm", "--phasecenter-remove", dest="phase_center_rm",
        type=float, default=None, help="The radius from the phase center"
        "not to consider for final direction-dependent source selection."
        "In beam sizes. Default is None.")

    add('-jc', '--json-config', dest='config', default=None,
        help='Json config file : No default')

    args = parser.parse_args()
    
    pybdsm_opts = dict([ items.split("=") for items in args.to_pybdsm ] ) \
                       if args.to_pybdsm else {}

    if args.rel_src_excl:
        rel_rmsrc = args.rel_src_excl[0].split(':')
    else:
        rel_rmsrc = None


    # making outdirectory
    def get_prefix(prefix, imagename, outdir):
        if prefix is None:
            prefix = os.path.basename(imagename.split(",")[0]).split(".")[:-1]
            prefix = ".".join(prefix)

        outdir = outdir or prefix +"_"+ os.path.join(datetime.datetime.now().\
                                           strftime("%Y-%m-%d-%H"))
        outdir = os.path.abspath(outdir)

        if not os.path.exists(outdir): 
            os.makedirs(outdir)

        prefix = outdir +"/"+ prefix
        
        return prefix
    
    if args.config:

        keys = []

        images = args.image.split(",") if args.image else None

        with open(args.config) as conf:
            jdict = json.load(conf)
            for key,val in jdict.items():
                if isinstance(val, unicode):
                    jdict[key] = str(val)
        
        prefix = get_prefix(jdict["prefix"],
                    images[0] if images else jdict["imagename"],
                    jdict["outdir"])

        reldict = jdict["reliability"]
        ddict = jdict["dd_tagging"]
        sourcefin = jdict["source_finder_opts"]

        for key, val in reldict.iteritems():
            if isinstance(val, unicode):
                reldict[key] = str(val)
        
        for key, val in ddict.iteritems():
            if isinstance(val, unicode):
                ddict[key] = str(val)


        reldict.update(sourcefin)

        enable = ddict["enable"]
        ddict.pop("enable")

        if args.image:
            psfs = args.psf.split(",") if args.psf else [None]
        
            if len(images) != len(psfs):
                psfs = [psfs[0]]*len(images)

            for i, (image, psf) in enumerate(zip(images, psfs)):
                if len(images) >1:
                    prefix = prefix + "-%04d"%i

                reldict["prefix"]  = prefix
                mc = rel.load(image, psf, **reldict)
                pos, neg = mc.get_reliability()

                ddict["poscatalog"] = pos
                ddict["negcatalog"] = neg
                ddict["prefix"] = prefix

                if enable and args.psf:
                    dc = dd.load(image, psf, **ddict)
                    ppos, nneg = dc.source_selection()

        else:

            image = jdict["imagename"]
            psf = jdict["psfname"]

            reldict["prefix"]  = prefix 
            mc = rel.load(image, psf, **reldict)
            pos, neg = mc.get_reliability()

            ddict["poscatalog"] = pos
            ddict["negcatalog"] = neg
            ddict["prefix"] = prefix

            if enable and psf:
                dc = dd.load(image, psf, **ddict)
                ppos, nneg = dc.source_selection()
        
    else:
        # reliability
        images = args.image.split(",")
        psfs = args.psf.split(",") if args.psf else [None]


        psfregion = args.psfcorr_region
        locregion = args.local_region

        if len(images) != len(psfs):
            psfs = [psfs[0]]*len(images)

        for i, (image, psf) in enumerate(zip(images, psfs)):
            prefix = get_prefix(args.prefix, images[0], args.outdir)

            if len(images) >1:
                prefix = prefix + "-%04d"%i

            mc  = rel.load(imagename=image, psfname=psf, sourcefinder_name=
                     args.source_finder, makeplots=args.do_relplots, 
                     do_psf_corr=args.add_psfcorr, do_local_var=args.add_locvar,
                     psf_corr_region=psfregion, local_var_region=locregion, 
                     rel_excl_src=rel_rmsrc, pos_smooth=args.pos_smooth,
                     neg_smooth=args.neg_smooth, loglevel=args.log_level, 
                     thresh_isl=args.thresh_isl, thresh_pix=args.thresh_pix,
                     neg_thresh_isl=args.neg_thresh_isl, neg_thresh_pix=
                     args.neg_thresh_pix, prefix=prefix, **pybdsm_opts)

            # assignign reliability values
            pos, neg = mc.get_reliability()

            # direction dependent detection tagging
            if args.psf:
                dc = dd.load(imagename=image, psfname=psf, poscatalog=pos, negcatalog=neg,
                           snr_thresh=args.snr_thresh, local_thresh=args.locvar_thresh,
                           local_region=locregion, psfcorr_region=psfregion, 
                           high_corr_thresh=args.psfcorr_thresh, negdetec_region=
                           args.neg_region, negatives_thresh=args.num_negatives,
                           phasecenter_excl_radius=args.phase_center_rm, prefix=prefix,
                           loglevel=args.log_level)
            # tagging
                ppose, nneg = dc.source_selection()

    os.system("rm -r *.pybdsm.log")
