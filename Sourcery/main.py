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
        help="Point Spread Function (PSF) Fits image name. Default = None")

    add("-s", "--source-finder", dest="source_finder", default="pybdsm",
        help="Source finder name. Default='pybdsm'")

    add("--to-pybdsm", dest="to_pybdsm", action="append",
           help="PyBDSM process_image options" 
           " [may be specified multiple times]. E.g thresh_pix=1 etc.")

    add("-log", "--log-level",dest="log_level", type=int, default=0,
        help="Python logging module, where info=0, debug=1, critial=2 and abort=3")

    add("-od", "--output-dir", dest="outdir", default=None, 
        help="Output products will be dumped here. System Default will be generated")

    add("-pref", "-dir-prefix", dest="prefix", default=None, 
        help="Prefix for output products. System default will be generated")

    add("-apsf", "--add-psfcorr", dest="add_psfcorr", action="store_true",
        default=False, help="Adds PSF correlation for density estimations. The PSF "
        " image must be provided. Default is False. To set True add -apsf")

    add("-alv", "--add-localvar", dest="add_locvar", action="store_true",
        default=False, help=" Include local variance for density estimations. See -apsf.")

    add("-dn", "--do-nearsources", dest="do_nearsources", action="store_true",
        default=False, help=" Adds number of neighbours for density estimations. See -apsf")

    add("-dmp", "--do-relplots", dest="do_relplots", action="store_false",
        default=True, help="Make reliability density plot. Default is True."
        " To disable add -dmp on the command line.")
 
    add("-rel", "--rel-thresh", dest="rel_thresh", default=None, type=float, 
        help= "Sets a reliability threshold. Default is None.")

    add("-beam", "--beam-cluster", dest="do_beam", default=False,
        action="store_true", help= "Increases the Gaussian groupings by 20 percent the"
        " beam size. Default is False.")
 
    add("-pcr", "--psfcorr-region",  dest="psfcorr_region", type=int,
        default=5, help="The size of the region to correlate, in beam sizes."
        " Default value is 5 beam sizes.")
  
    add("-lr", "--local-region", dest="local_region", default=10,
        help="The area to compute the local variance, in beam sizes."
        " Default value is 10 beam sizes.")
     
    add("-rel_rm", "--rel-sources-excl", dest="rel_src_excl", 
        action="append", default=None, help="Delete sources within a radius;"
        " e.g ra, dec, radius (in degrees). For more than"
        " one region: ra1,dec1,radius1:ra2,dec2,radius2. Default is None.")

    add("-ps", "--positive-smooth", dest="pos_smooth", type=float, default=2.0, 
        help=" Masking threshold in the positive image, e.g pixels 2x noise"
        "will be masked. Default is 2.0.")

    add("-ns", "--negative-smooth", dest="neg_smooth", type=float, default=2.0, 
        help="Similar to -ps. Applied to negative image")

    add('-pisl', "--thresh-isl", dest="thresh_isl", type=float, default=3,
        help=" Source finding threshold for forming islands. For positive image."
         "  Default is 3.")

    add("-ppix", "--thresh-pix", dest="thresh_pix", type=float, default=5,
        help="Source finding threshold for model fitting. For positive image."
        "Default is 5.")
 
    add("-nisl", "--negthreshold-isl", dest="neg_thresh_isl",
        type=float, default=3, help="Similar to -pisl but applied"
        " to the negative image. Default is 3.")

    add("-npix", "--negthreshold-pix", dest="neg_thresh_pix",
        type=float, default=5, help="Similar to -ppix but applied"
        " for negative image. Default is 5.")

    add("-snr_thr", "--snr-threshold", dest="snr_thresh", type=float,
        default=40, help="SNR threshold. High SNR Sources have SNR > 40x min(SNR)"
        " Default is 40.")

    add("-loc_thr", "--localvar-threshold", dest="locvar_thresh",
        type=float, default=0.4, help="Local variance (LV) threshold."
        "High LV Sources have LV > 0.4 x image noise. Default is 0.4")

    add("-pc_thr", "--psfcorr-threshold", dest="psfcorr_thresh", 
        type=float, default=0.5, help="Correlation factor (CF) threshold."
        "High CF sources have CF > 0.5. Default is 0.5.")    

    add("-nneg", "--num-negatives",dest="num_negatives", type=float, default=4,
        help="Number of neigbouring (NN) negative detections. Default is 4."
        "High NN sources have NN > 4.")

    add("-nrgn", "--neg-region", dest="neg_region", type=float, default=10,
        help="The size of a region to find -nneg. Default is 10.") 

    add("-nphrm", "--phasecenter-remove", dest="phase_center_rm",
        type=float, default=None, help="The radius excluded from"
        " direction-dependent source selection. NB: this radius is wrt to"
        " the phase center. Default is None.")

    add('-jc', '--json-config', dest='config', default=None,
        help='Json config file : No default')
    
    add("-smp", "--save-posmask", dest="savemask_pos", action="store_true",
        default=False, help="If specified, positive mask is saved.")
    
    add("-smn", "--save-negmask", dest="savemask_neg", action="store_true",
        default=False, help=" If specified, the negative mask is saved.")


    args = parser.parse_args() 
    pybdsm_opts = dict([ items.split("=") for items in args.to_pybdsm ] ) \
                       if args.to_pybdsm else {}

    if args.rel_src_excl:
        rel_rmsrc = args.rel_src_excl[0].split(':')
    else:
        rel_rmsrc = None

    if not args.image:
        print("ATTENTION: No image provided. Aborting")

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
                prefix = get_prefix(jdict["prefix"],
                         images[0] if images else jdict["imagename"],
                         jdict["outdir"])
                if len(images) >1:
                    prefix = prefix + "-%04d"%i

                reldict["prefix"]  = prefix
                mc = rel.load(image, psf, **reldict)
                pmodel, nmodel, noise, hdr = mc.get_reliability()

                ddict["pmodel"] = pmodel
                ddict["nmodel"] = nmodel
                ddict["prefix"] = prefix
                ddict["noise"] = noise
                ddict["header"] = hdr

                if enable:
                    dc = dd.load(image, psfname=psf, **ddict)
                    pmodel, nmodel  = dc.source_selection()

        else:

            image = jdict["imagename"]
            psf = jdict["psfname"]

            reldict["prefix"]  = prefix 
            mc = rel.load(image, psf, **reldict)
            pmodel, nmodel, noise, bmaj = mc.get_reliability()

            ddict["pmodel"] = pmodel
            ddict["nmodel"] = nmodel
            ddict["prefix"] = prefix
            ddict["noise"]  = noise
            ddict["header"] = hdr

            if enable:
                dc = dd.load(image, psfname=psf, **ddict)
                pmodel, nmodel = dc.source_selection()

        pmodel.save( prefix+".lsm.html")
        nmodel.save( prefix+"_negative.lsm.html")
        
    else:
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
                     args.neg_thresh_pix, prefix=prefix,  
                     do_nearsources=args.do_nearsources, increase_beam_cluster=
                     args.do_beam, savemask_neg=args.savemask_neg,
                     savemask_pos=args.savemask_pos, **pybdsm_opts)

            # assignign reliability values
            pmodel, nmodel, noise, hdr = mc.get_reliability()

            # direction dependent detection tagging
            
            dc = dd.load(imagename=image, psfname=psf, pmodel=pmodel, nmodel=nmodel,
                    header=hdr, snr_thresh=args.snr_thresh, local_thresh=args.locvar_thresh, 
                    high_corr_thresh=args.psfcorr_thresh, negdetec_region=
                    args.neg_region, negatives_thresh=args.num_negatives, noise=noise,
                    phasecenter_excl_radius=args.phase_center_rm, prefix=prefix,
                    loglevel=args.log_level)
            # tagging
            pmodel, nmodel = dc.source_selection()
            pmodel.save( prefix+".lsm.html")
            nmodel.save( prefix+"_negative.lsm.html")

    os.system("rm -r tmp*.log")
