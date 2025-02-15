#!/bin/sh





python ./scripts/run_tis_80310229001.py >  tis_out_80310229001.txt &
python ./scripts/run_tis_80310230001.py >  tis_out_80310230001.txt &
wait

echo "all orbit scripts finished"