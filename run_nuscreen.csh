#!/bin/csh
# usage: run_nuscreen.csh PATH/TO/OBSID
# Variables below must correspond to names of files

set interval=22-22-00_22-27-00
set fpm=A
set INDIR=/Users/jessieduncan/nustar/may-2018/5_29pixpos/80410205001/
set work_dir=./initial_dem/

# syntax: run_screensplit.sh INDIR
# if [ $# != 1 ] ; then
#     echo Syntax:  screensplitproducts.sh PATH_TO_OBSID
#     exit 1
# fi

set type="STATUS==b0000xx00xx0xx000"


set STEMINPUTS=nu`basename ${INDIR}`

nuscreen obsmode=SCIENCE_SC infile=$INDIR/event_cl/$STEMINPUTS$fpm"_uf.evt" \
gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT \
gradeexpr=0 statusexpr=$type \
createattgti=yes createinstrgti=yes \
clobber=yes \
outdir=$work_dir/$interval/ usrgtifile=$work_dir/$interval/$interval$fpm"_gti.fits" \
hkfile=$INDIR/hk/$STEMINPUTS$fpm"_fpm.hk" mkffile=$INDIR/event_cl/$STEMINPUTS$fpm".mkf" outfile=DEFAULT

mv $work_dir/$interval/$STEMINPUTS$fpm"06_cl.evt" $work_dir/$interval/$STEMINPUTS$fpm"06_0_p_cl.evt"

nuscreen obsmode=SCIENCE_SC infile=$INDIR/event_cl/$STEMINPUTS$fpm"_uf.evt" \
gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT \
gradeexpr=21-24 statusexpr=$type \
createattgti=yes createinstrgti=yes \
clobber=yes \
outdir=$work_dir/$interval/ usrgtifile=$work_dir/$interval/$interval$fpm"_gti.fits" \
hkfile=$INDIR/hk/$STEMINPUTS$fpm"_fpm.hk" mkffile=$INDIR/event_cl/$STEMINPUTS$fpm".mkf" outfile=DEFAULT

mv $work_dir/$interval/$STEMINPUTS$fpm"06_cl.evt" $work_dir/$interval/$STEMINPUTS$fpm"06_21_24_p_cl.evt"

nuscreen obsmode=SCIENCE_SC infile=$INDIR/event_cl/$STEMINPUTS$fpm"_uf.evt" \
gtiscreen=yes evtscreen=yes gtiexpr=DEFAULT \
gradeexpr=0-4 statusexpr=$type \
createattgti=yes createinstrgti=yes \
clobber=yes \
outdir=$work_dir/$interval/ usrgtifile=$work_dir/$interval/$interval$fpm"_gti.fits" \
hkfile=$INDIR/hk/$STEMINPUTS$fpm"_fpm.hk" mkffile=$INDIR/event_cl/$STEMINPUTS$fpm".mkf" outfile=DEFAULT

mv $work_dir/$interval/$STEMINPUTS$fpm"06_cl.evt" $work_dir/$interval/$STEMINPUTS$fpm"06_0_4_p_cl.evt"
