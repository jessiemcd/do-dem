#!/bin/sh





python ./scripts/run_tis_20617005001.py >  tis_out_20617005001.txt &
python ./scripts/run_tis_20617006001.py >  tis_out_20617006001.txt &
wait

echo "all orbit scripts finished"