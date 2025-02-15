#!/bin/sh





python ./scripts/run_tis_20611005001.py >  tis_out_20611005001.txt &
python ./scripts/run_tis_20611006001.py >  tis_out_20611006001.txt &
wait

echo "all orbit scripts finished"