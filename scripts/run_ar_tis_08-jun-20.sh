#!/bin/sh





python ./scripts/run_tis_20611008001.py >  tis_out_20611008001.txt &
python ./scripts/run_tis_20611009001.py >  tis_out_20611009001.txt &
wait

echo "all orbit scripts finished"