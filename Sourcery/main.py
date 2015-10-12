# Reliability estimation and direction-dependent source selection


import sys
from argparse import ArgumentParser
import reliabilityestimates as rel
import direcdepen as dd
import os
import datetime
import json

def main():

    __version_info__ = (0,0,1)
    __version__ = ".".join( map(str,__version_info__) )

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

    add("-pref", "-dir-prefix", dest="prefix", default=None, 
        help="Give a prefix to an output directory. Default is None")

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
        help="Data size to compute the local variance, in beam sizes."
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

    add('-jc', '--json-config', dest='config',
        help='Json config file : No default')

    args = parser.parse_args()
    
    # image and psf image
    img = args.image
    psf = args.psf

    #source finder
    sf = args.source_finder
    psmth = args.pos_smooth
    nsmth = args.neg_smooth
    nthresh_isl = args.neg_thresh_isl
    nthresh_pix = args.neg_thresh_pix
    pthr_isl = args.thresh_isl
    pthr_pix = args.thresh_pix
    pybdsm_opts = dict([ items.split("=") for items in args.to_pybdsm ] ) \
                       if args.to_pybdsm else {}

    #boolean options
    mkplt = args.do_relplots
    dopsf = args.add_psfcorr
    doloc = args.add_locvar
        
    # log level
    log = args.log_level

    # removing sources at some region
    if args.rel_src_excl:
        rel_rmsrc = args.rel_src_excl[0].split(':')
    else:
        rel_rmsrc = None


    # making outdirectory
    if args.prefix is None:
        prefix = os.path.basename(args.image).split(".")[:-1]
        prefix = ".".join(prefix)
    else:
        prefix = args.prefix

    outdir = prefix +"_"+ os.path.join(datetime.datetime.now().\
                                       strftime("%Y-%m-%d-%H"))
    outdir = os.path.abspath(outdir)

    if not os.path.exists(outdir): 
        os.makedirs(outdir)

    prefix = outdir +"/"+ prefix

    # thresholds
    sthr = args.snr_thresh
    lthr = args.locvar_thresh
    pcthr = args.psfcorr_thresh
    negthr = args.num_negatives

    # regions to evaluate
    psfregion = args.psfcorr_region
    locregion = args.local_region
    negregion = args.neg_region
    phase_remove = args.phase_center_rm
  
    
    if args.config:
        keys = []
        with open(args.config) as conf:
            jdict = json.load(conf)
        
        reldict = jdict["reliability"]
        ddict = jdict["dd_tagging"]
        sourcefin = jdict["source_finder_opts"]

        for key, val in reldict.iteritems():
            if isinstance(val, unicode):
                reldict[key] = str(val)
        
        z = reldict.copy()
        z.update({"prefix" : prefix})
        z.update(sourcefin)
        mc = rel.compute(**z)
        pos, neg = mc.get_reliability()

        for key, val in ddict.iteritems():
            if isinstance(val, unicode):
                ddict[key] = str(val)

        x = ddict.copy()
        y = {"poscatalog" : pos, "negcatalog" : neg, "prefix" : prefix}
        x.update(y)
        dc = dd.Parent(**x)
        ppos, nneg = dc.source_selection()
        
    else:
        # reliability     
        mc  = rel.compute(imagename=img, psfname=psf, sourcefinder_name=sf,
                 makeplots=mkplt, do_psf_corr=dopsf, do_local_var=doloc,
                 psf_corr_region=psfregion, local_var_region=locregion, 
                 rel_excl_src=rel_rmsrc, pos_smooth=psmth,
                 neg_smooth=nsmth, loglevel=log, thresh_isl=pthr_isl,
                 thresh_pix=pthr_pix, neg_thresh_isl=nthresh_isl,
                 neg_thresh_pix=nthresh_pix, prefix=prefix, **pybdsm_opts)

        # assignign reliability values
        pos, neg = mc.get_reliability()

        # direction dependent detection tagging
        dc = dd.Parent(imagename=img, psfname=psf, poscatalog=pos, negcatalog=neg,
                   snr_thresh=sthr, local_thresh=lthr, local_region=locregion, 
                   psfcorr_region=psfregion, high_corr_thresh=pcthr,
                   negdetec_region=negregion, negatives_thresh=negthr,
                   phasecenter_excl_radius=phase_remove, prefix=prefix,
                   loglevel=log)
        # tagging
        ppose, nneg = dc.source_selection()
