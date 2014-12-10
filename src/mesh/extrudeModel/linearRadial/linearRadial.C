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

#include "linearRadial.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linearRadial, 0);

addToRunTimeSelectionTable(extrudeModel, linearRadial, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearRadial::linearRadial(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    R_(readScalar(coeffDict_.lookup("R"))),
    Rsurface_(coeffDict_.lookupOrDefault<scalar>("Rsurface", -1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearRadial::~linearRadial()
{}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point linearRadial::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    // radius of the surface
    scalar rs = mag(surfacePoint);
    vector rsHat = surfacePoint/rs;
    if (Rsurface_ >= 0) rs = Rsurface_;

    scalar r = rs + (R_ - rs)*sumThickness(layer);
    return r*rsHat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //
