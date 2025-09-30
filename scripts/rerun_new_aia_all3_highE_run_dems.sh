#!/bin/sh



python ./scripts/dems_26-jul-16_1.py >  ./scripts/dems_out_26-jul-16_1.txt &
python ./scripts/dems_19-feb-16.py >  ./scripts/dems_out_19-feb-16.txt &
python ./scripts/dems_02-sep-15.py >  ./scripts/dems_out_02-sep-15.txt &
wait

echo "all orbit scripts finished"