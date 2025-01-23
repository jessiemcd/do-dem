#!/bin/sh
# usage: run_nuproducts.sh
# Variables below must correspond to names of files

interval=14-59-45_15-03-20
region=/Users/jmdunca2/do-dem/initial_dem_apr21_gauss/14-59-45_15-03-20/nu20615001001A06_0_4_p_cl_sunpos__region.reg
fpm=A
INDIR=/Users/jmdunca2/nustar/apr-2021/20615001001/
unphys_products=1
working_dir=/Users/jmdunca2/do-dem/initial_dem_apr21_gauss/
adjacent_grades=1

type="STATUS==b0000xx00xx0xx000"


STEMINPUTS=nu`basename ${INDIR}`


if [ $adjacent_grades == 0 ] ; then
    echo Making spectral data products for grade 0

    nuproducts indir=$INDIR/event_cl/ \
    instrument=FPM$fpm \
    steminputs=$STEMINPUTS \
    outdir=$working_dir/$interval \
    extended=no \
    runmkarf=yes \
    runmkrmf=yes \
    infile=$working_dir/$interval/$STEMINPUTS$fpm"06_0_p_cl.evt" \
    bkgextract=no \
    srcregionfile=$region \
    attfile=$INDIR/event_cl/$STEMINPUTS"_att.fits" \
    hkfile=$INDIR/event_cl/$STEMINPUTS$fpm"_fpm.hk" \
    clobber=yes

fi


if [ $adjacent_grades == 1 ] ; then
    echo Making spectral data products for grades 0-4

    nuproducts indir=$INDIR/event_cl/ \
    instrument=FPM$fpm \
    steminputs=$STEMINPUTS \
    outdir=$working_dir/$interval \
    extended=no \
    runmkarf=yes \
    runmkrmf=yes \
    infile=$working_dir/$interval/$STEMINPUTS$fpm"06_0_4_p_cl.evt" \
    bkgextract=no \
    srcregionfile=$region \
    attfile=$INDIR/event_cl/$STEMINPUTS"_att.fits" \
    hkfile=$INDIR/event_cl/$STEMINPUTS$fpm"_fpm.hk" \
    clobber=yes

fi

if [ $unphys_products != 1 ] ; then
    echo Not making spectral products for unphysical grades
    exit 1
fi

nuproducts indir=$INDIR/event_cl/ \
instrument=FPM$fpm \
steminputs=$STEMINPUTS \
outdir=$working_dir/$interval \
extended=no \
runmkarf=yes \
runmkrmf=yes \
infile=$working_dir/$interval/$STEMINPUTS$fpm"06_21_24_p_cl.evt" \
bkgextract=no \
srcregionfile=$region \
attfile=$INDIR/event_cl/$STEMINPUTS"_att.fits" \
hkfile=$INDIR/event_cl/$STEMINPUTS$fpm"_fpm.hk" \
clobber=yes


