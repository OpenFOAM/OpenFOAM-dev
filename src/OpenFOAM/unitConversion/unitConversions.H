/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

Global
    unitConversions

Description
    Useful unit conversions

SourceFiles
    unitConversions.C

\*---------------------------------------------------------------------------*/

#ifndef unitConversions_H
#define unitConversions_H

#include "unitConversion.H"
#include "dimensionSets.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern const unitConversion unitless;

extern const unitConversion unitAny;
extern const unitConversion unitNone;

extern const unitConversion unitFraction;
extern const unitConversion unitPercent;

extern const unitConversion unitRadians;
extern const unitConversion unitRotations;
extern const unitConversion unitDegrees;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Convert degrees to radians
scalar degToRad(const scalar deg);

//- Convert radians to degrees
scalar radToDeg(const scalar rad);

//- Add a unit conversion
void addUnits(const word& name, const unitConversion& units);

//- Get the table of unit conversions
const HashTable<unitConversion>& units();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
