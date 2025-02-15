#!/bin/sh





python ./scripts/run_tis_20201001001.py >  tis_out_20201001001.txt &
python ./scripts/run_tis_20201002001.py >  tis_out_20201002001.txt &
python ./scripts/run_tis_20201003001.py >  tis_out_20201003001.txt &
python ./scripts/run_tis_20201004001.py >  tis_out_20201004001.txt &
wait

echo "all orbit scripts finished"