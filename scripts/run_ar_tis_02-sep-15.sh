#!/bin/sh





python ./scripts/run_tis_20102005001.py >  tis_out_20102005001.txt &
python ./scripts/run_tis_20102006001.py >  tis_out_20102006001.txt &
python ./scripts/run_tis_20102007001.py >  tis_out_20102007001.txt &
python ./scripts/run_tis_20102008001.py >  tis_out_20102008001.txt &
wait

echo "all orbit scripts finished"