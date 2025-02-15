#!/bin/sh





python ./scripts/run_tis_80416208001.py >  tis_out_80416208001.txt &
python ./scripts/run_tis_80416209001.py >  tis_out_80416209001.txt &
python ./scripts/run_tis_80416210001.py >  tis_out_80416210001.txt &
python ./scripts/run_tis_80416211001.py >  tis_out_80416211001.txt &
python ./scripts/run_tis_80416212001.py >  tis_out_80416212001.txt &
python ./scripts/run_tis_80416213001.py >  tis_out_80416213001.txt &
wait

echo "all orbit scripts finished"