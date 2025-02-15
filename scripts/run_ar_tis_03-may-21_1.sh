#!/bin/sh





python ./scripts/run_tis_20616001001.py >  tis_out_20616001001.txt &
python ./scripts/run_tis_20616002001.py >  tis_out_20616002001.txt &
wait

echo "all orbit scripts finished"