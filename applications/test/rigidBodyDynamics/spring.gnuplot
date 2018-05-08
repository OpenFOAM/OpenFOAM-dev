#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
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
#     spring.gnuplot
#
# Description
#     Creates an PostScript graph file of Test-spring results vs
#     the analytical solution.
#
#------------------------------------------------------------------------------

reset

set samples 2000

k = 5000.0
m = 9.6
c = 50.0
a = -0.1

omega = sqrt(k/m)
zeta = c/(2.0*m*omega)

phi = atan((sqrt(1.0 - zeta**2))/zeta)
A = a/sin(phi)

pos(A, t, omega, phi, zeta) = A*exp(-zeta*omega*t)*sin(sqrt(1-zeta**2)*omega*t + phi)
vel(A, t, omega, phi, zeta) = \
A*exp(-zeta*omega*t)*\
( \
  sqrt(1-zeta**2)*omega*cos(sqrt(1-zeta**2)*omega*t + phi) \
- zeta*omega*sin(sqrt(1-zeta**2)*omega*t + phi) \
)

set xlabel "Time (s)"
set ylabel "Position (m)"
set y2label "Velocity (m/s)"

set ytics nomirror
set y2tics

set yrange [-0.1:0.1]
set y2range [-2:2]

set xzeroaxis

set terminal postscript eps color enhanced solid
set output "spring.eps"

plot \
    "spring.dat" u 1:($2 - 0.1) w l t "Simulation, centre of mass relative to start", \
    pos(A, x, omega, phi, zeta) w l t "Analytical solution, centre of mass", \
    "spring.dat" u 1:3 w l axes x1y2 t "Simulation, vertical velocity", \
    vel(A, x, omega, phi, zeta) w l axes x1y2 t "Analytical solution, vertical velocity"

#------------------------------------------------------------------------------
