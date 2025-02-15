#!/bin/sh





python ./scripts/run_tis_80416201001.py >  tis_out_80416201001.txt &
python ./scripts/run_tis_80416202001.py >  tis_out_80416202001.txt &
python ./scripts/run_tis_80416203001.py >  tis_out_80416203001.txt &
python ./scripts/run_tis_80416204001.py >  tis_out_80416204001.txt &
python ./scripts/run_tis_80416205001.py >  tis_out_80416205001.txt &
wait

echo "all orbit scripts finished"