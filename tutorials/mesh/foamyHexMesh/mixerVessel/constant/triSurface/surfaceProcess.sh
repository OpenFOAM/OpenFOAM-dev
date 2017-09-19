#!/bin/sh

cd ${0%/*} || exit 1    # Run from this directory

cp rawSurfaces/*.stl .


# Vessel surface
surfaceAdd outlet.stl vessel.stl vessel.stl

# Sparger
surfaceCheck sparger.stl
surfaceAdd gasInlet.stl sparger_0.obj spargerInlet.stl
surfaceConvert sparger_1.obj spargerShaft.stl
surfaceOrient -inside spargerInlet.stl spargerInlet.stl "(0 0.1 1)"
surfaceOrient -inside spargerShaft.stl spargerShaft.stl "(0 0.1 -0.1)"

# Rotating shaft
surfaceOrient -inside shaftRotating.stl shaftRotating.stl "(0 0.1 1)"

# Static shaft
surfaceOrient -inside shaft.stl shaft.stl "(0 0.1 1)"

# Stirrer
surfaceSplitByTopology stirrer.stl stirrer.stl
surfaceOrient -inside stirrer.stl stirrer.stl "(0 0.1 1)"

mv stirrer_bafflePart_0.stl stirrer_baffles.stl

surfaceCheck stirrer_baffles.stl
mv stirrer_baffles_0.obj stirrer_baffles_plate.obj


#------------------------------------------------------------------------------
