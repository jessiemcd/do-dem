#!/bin/sh





python ./scripts/run_tis_20101017001.py >  tis_out_20101017001.txt &
python ./scripts/run_tis_20101018001.py >  tis_out_20101018001.txt &
python ./scripts/run_tis_20101019001.py >  tis_out_20101019001.txt &
python ./scripts/run_tis_20101020001.py >  tis_out_20101020001.txt &
wait

echo "all orbit scripts finished"