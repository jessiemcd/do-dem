#!/bin/sh





python ./scripts/run_tis_20615001001.py >  tis_out_20615001001.txt &
python ./scripts/run_tis_20615002001.py >  tis_out_20615002001.txt &
python ./scripts/run_tis_20615003001.py >  tis_out_20615003001.txt &
python ./scripts/run_tis_20615004001.py >  tis_out_20615004001.txt &
python ./scripts/run_tis_20615005001.py >  tis_out_20615005001.txt &
wait

echo "all orbit scripts finished"