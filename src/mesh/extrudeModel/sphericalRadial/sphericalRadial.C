/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "sphericalRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sphericalRadial, 0);

addToRunTimeSelectionTable(extrudeModel, sphericalRadial, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sphericalRadial::sphericalRadial(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    R_(Function1<scalar>::New("R", coeffDict_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

sphericalRadial::~sphericalRadial()
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point sphericalRadial::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    // Radius of the surface
    const scalar rs = mag(surfacePoint);

    // Radial direction of surfacePoint
    const vector rsHat = surfacePoint/rs;

    // Radius of layer
    const scalar r = R_->value(layer);

    // Return new point
    return r*rsHat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //
