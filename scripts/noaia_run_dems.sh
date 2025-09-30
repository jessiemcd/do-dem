#!/bin/sh



python ./scripts/dems_08-jun-20.py >  ./scripts/dems_out_08-jun-20.txt &
python ./scripts/dems_08-jan-21.py >  ./scripts/dems_out_08-jan-21.txt &
python ./scripts/dems_26-jul-16_2.py >  ./scripts/dems_out_26-jul-16_2.txt &
wait

echo "all orbit scripts finished"