#!/bin/sh





python ./scripts/run_tis_20101015001.py >  tis_out_20101015001.txt &
python ./scripts/run_tis_20101016001.py >  tis_out_20101016001.txt &
wait

echo "all orbit scripts finished"