#!/bin/sh





python ./scripts/run_tis_20102002001.py >  tis_out_20102002001.txt &
python ./scripts/run_tis_20102003001.py >  tis_out_20102003001.txt &
python ./scripts/run_tis_20102004001.py >  tis_out_20102004001.txt &
wait

echo "all orbit scripts finished"