#!/bin/sh





python ./scripts/run_tis_80310211001.py >  tis_out_80310211001.txt &
python ./scripts/run_tis_80310212001.py >  tis_out_80310212001.txt &
python ./scripts/run_tis_80310213001.py >  tis_out_80310213001.txt &
python ./scripts/run_tis_80310214001.py >  tis_out_80310214001.txt &
wait

echo "all orbit scripts finished"