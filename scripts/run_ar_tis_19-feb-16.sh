#!/bin/sh





python ./scripts/run_tis_20102011001.py >  tis_out_20102011001.txt &
python ./scripts/run_tis_20102012001.py >  tis_out_20102012001.txt &
python ./scripts/run_tis_20102013001.py >  tis_out_20102013001.txt &
python ./scripts/run_tis_20102014001.py >  tis_out_20102014001.txt &
wait

echo "all orbit scripts finished"