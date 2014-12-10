/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "wedge.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(wedge, 0);

addToRunTimeSelectionTable(extrudeModel, wedge, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wedge::wedge(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    axisPt_(coeffDict_.lookup("axisPt")),
    axis_(coeffDict_.lookup("axis")),
    angle_
    (
        degToRad(readScalar(coeffDict_.lookup("angle")))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

wedge::~wedge()
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point wedge::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    scalar sliceAngle;
    // For the case of a single layer extrusion assume a
    // symmetric wedge about the reference plane is required
    if (nLayers_ == 1)
    {
        if (layer == 0)
        {
            sliceAngle = -angle_/2.0;
        }
        else
        {
            sliceAngle = angle_/2.0;
        }
    }
    else
    {
        //sliceAngle = angle_*layer/nLayers_;
        sliceAngle = angle_*sumThickness(layer);
    }

    // Find projection onto axis (or rather decompose surfacePoint
    // into vector along edge (proj), vector normal to edge in plane
    // of surface point and surface normal.
    point d = surfacePoint - axisPt_;

    d -= (axis_ & d)*axis_;

    scalar dMag = mag(d);

    point edgePt = surfacePoint - d;

    // Rotate point around sliceAngle.
    point rotatedPoint = edgePt;

    if (dMag > VSMALL)
    {
        vector n = (d/dMag) ^ axis_;

        rotatedPoint +=
          + cos(sliceAngle)*d
          - sin(sliceAngle)*mag(d)*n; // Use either n or surfaceNormal
    }

    return rotatedPoint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //
