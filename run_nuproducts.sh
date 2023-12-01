#!/bin/sh
# usage: run_nuproducts.sh
# Variables below must correspond to names of files

interval=22-22-00_22-32-00
region=initial_dem/22-22-00_22-32-00/nu80410205001B06_0_p_cl_sunpos_COM_region.reg
fpm=B
INDIR=/Users/jessieduncan/nustar/may-2018/5_29pixpos/80410205001/
unphys_products=0
working_dir=./initial_dem/

type="STATUS==b0000xx00xx0xx000"


STEMINPUTS=nu`basename ${INDIR}`


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


