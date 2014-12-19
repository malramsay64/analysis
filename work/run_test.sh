#!/bin/sh

#  Script.sh
#  analysis
#
#  Created by Malcolm Ramsay on 19/12/2014.
#  Copyright (c) 2014 Malcolm Ramsay. All rights reserved.

gnuplot=/usr/local/bin/gnuplot

#echo hello world
$gnuplot unit_cell.plot
#gnuplot -e "filename='trj_contact/0000000000'" ../gnuplot/frame.plot