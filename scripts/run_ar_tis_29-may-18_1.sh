#!/bin/sh





python ./scripts/run_tis_80410201001.py >  tis_out_80410201001.txt &
python ./scripts/run_tis_80410202001.py >  tis_out_80410202001.txt &
python ./scripts/run_tis_80410203001.py >  tis_out_80410203001.txt &
python ./scripts/run_tis_80410204001.py >  tis_out_80410204001.txt &
python ./scripts/run_tis_80410205001.py >  tis_out_80410205001.txt &
wait

echo "all orbit scripts finished"