#!/bin/sh





python ./scripts/run_tis_20513002001.py >  tis_out_20513002001.txt &
python ./scripts/run_tis_20513005001.py >  tis_out_20513005001.txt &
python ./scripts/run_tis_20513006001.py >  tis_out_20513006001.txt &
python ./scripts/run_tis_20513007001.py >  tis_out_20513007001.txt &
python ./scripts/run_tis_20513008001.py >  tis_out_20513008001.txt &
wait

echo "all orbit scripts finished"