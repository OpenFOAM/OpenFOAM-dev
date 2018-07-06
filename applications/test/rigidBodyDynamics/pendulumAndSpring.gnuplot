#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     | Website:  https://openfoam.org
#   \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
#    \\/     M anipulation  |
#-------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     pendulumAndSpring.gnuplot
#
# Description
#     Creates an PostScript graph file of Test-pendulumAndSpring results
#
#------------------------------------------------------------------------------

reset

set xlabel "Time/[s]"
set ylabel "x"
set y2label "omega"

set ytics nomirror
set y2tics

set yrange [-1.5:1.5]
set y2range [-35:35]

set xzeroaxis

set terminal postscript eps color enhanced solid
set output "pendulumAndSpring.eps"

plot \
    "pendulumAndSpring.dat" u 1:2 w l t "x", \
    "pendulumAndSpring.dat" u 1:(57.29578*$2) w l axes x1y2 t "omega"

#------------------------------------------------------------------------------
