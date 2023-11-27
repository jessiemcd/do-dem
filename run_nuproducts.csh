#!/bin/csh
# usage: run_nuproducts.csh PATH/TO/OBSID
# Variables below must correspond to names of files

set interval=22-22-00_22-27-00
set region=initial_dem/22-22-00_22-27-00/nu80410205001A06_0_p_cl_sunpos_COM_region.reg
set fpm=A
set INDIR=/Users/jessieduncan/nustar/may-2018/5_29pixpos/80410205001/
set unphys_products=1
set working_dir=./initial_dem/

set type="STATUS==b0000xx00xx0xx000"


set STEMINPUTS=nu`basename ${INDIR}`


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


