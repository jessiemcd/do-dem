#!/bin/sh
# usage: run_nuscreen.sh
# Variables below must correspond to names of files

interval=18-41-20_18-42-25
fpm=B
INDIR=/Users/jmdunca2/nustar/apr-2021/20615003001/
work_dir=/Users/jmdunca2/do-dem/initial_dem_apr21_gauss/
unphys_products=1
adjacent_grades=1

# syntax: run_screensplit.sh INDIR
# if [ $# != 1 ] ; then
#     echo Syntax:  screensplitproducts.sh PATH_TO_OBSID
#     exit 1
# fi

type="STATUS==b0000xx00xx0xx000"


STEMINPUTS=nu`basename ${INDIR}`

if [ $adjacent_grades == 0 ] ; then

    nuscreen obsmode=SCIENCE_SC infile=$INDIR/event_cl/$STEMINPUTS$fpm"_uf.evt" \
    gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT \
    gradeexpr=0 statusexpr=$type \
    createattgti=yes createinstrgti=yes \
    clobber=yes \
    outdir=$work_dir/$interval/ usrgtifile=$work_dir/$interval/$interval$fpm"_gti.fits" \
    hkfile=$INDIR/hk/$STEMINPUTS$fpm"_fpm.hk" mkffile=$INDIR/event_cl/$STEMINPUTS$fpm".mkf" outfile=DEFAULT

    mv $work_dir/$interval/$STEMINPUTS$fpm"06_cl.evt" $work_dir/$interval/$STEMINPUTS$fpm"06_0_p_cl.evt"

fi

if [ $unphys_products == 1 ] ; then

    nuscreen obsmode=SCIENCE_SC infile=$INDIR/event_cl/$STEMINPUTS$fpm"_uf.evt" \
    gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT \
    gradeexpr=21-24 statusexpr=$type \
    createattgti=yes createinstrgti=yes \
    clobber=yes \
    outdir=$work_dir/$interval/ usrgtifile=$work_dir/$interval/$interval$fpm"_gti.fits" \
    hkfile=$INDIR/hk/$STEMINPUTS$fpm"_fpm.hk" mkffile=$INDIR/event_cl/$STEMINPUTS$fpm".mkf" outfile=DEFAULT
    
    mv $work_dir/$interval/$STEMINPUTS$fpm"06_cl.evt" $work_dir/$interval/$STEMINPUTS$fpm"06_21_24_p_cl.evt"

fi

if [ $adjacent_grades == 1 ] ; then

    nuscreen obsmode=SCIENCE_SC infile=$INDIR/event_cl/$STEMINPUTS$fpm"_uf.evt" \
    gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT \
    gradeexpr=0-4 statusexpr=$type \
    createattgti=yes createinstrgti=yes \
    clobber=yes \
    outdir=$work_dir/$interval/ usrgtifile=$work_dir/$interval/$interval$fpm"_gti.fits" \
    hkfile=$INDIR/hk/$STEMINPUTS$fpm"_fpm.hk" mkffile=$INDIR/event_cl/$STEMINPUTS$fpm".mkf" outfile=DEFAULT
    
    mv $work_dir/$interval/$STEMINPUTS$fpm"06_cl.evt" $work_dir/$interval/$STEMINPUTS$fpm"06_0_4_p_cl.evt"

fi