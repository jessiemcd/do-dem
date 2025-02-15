#!/bin/sh





python ./scripts/run_tis_20611002001.py >  tis_out_20611002001.txt &
python ./scripts/run_tis_20611003001.py >  tis_out_20611003001.txt &
wait

echo "all orbit scripts finished"