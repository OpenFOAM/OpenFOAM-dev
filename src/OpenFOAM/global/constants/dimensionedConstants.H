/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    dimensionedConstants

Description
    Dictionary reading and supplying the dimensioned constants used within
    OpenFOAM, particularly for thermodynamics.

    The values are read from the OpenFOAM etc/controlDict and should be
    changed to run with a different set of units from the default SI units.

SourceFiles
    dimensionedConstants.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedConstants_H
#define dimensionedConstants_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Construct and return a constant by looking up from the DimensionedConstants
//  sub-dictionary of controlDict
dimensionedScalar dimensionedConstant
(
    const char* const group,
    const char* name,
    const dimensionSet& dimensions
);

//- Construct a dimensioned constant with a given value and insert into the
//  DimensionedConstants sub-dictionary of controlDict
dimensionedScalar dimensionedConstant
(
    const char* const group,
    const char* name,
    const dimensionedScalar& value
);

//- Construct a dimensioned constant with a given value and insert into the
//  DimensionedConstants sub-dictionary of controlDict with a different name
dimensionedScalar dimensionedConstant
(
    const char* const group,
    const char* entryName,
    const char* codeName,
    const dimensionedScalar& value
);

//- Construct a dimensioned constant with given units and value and insert into
//  the DimensionedConstants sub-dictionary of controlDict
dimensionedScalar dimensionedConstant
(
    const char* const group,
    const char* name,
    const unitConversion& units,
    const scalar value
);

//- Construct a dimensioned constant with given units and value and insert into
//  the DimensionedConstants sub-dictionary of controlDict with a different name
dimensionedScalar dimensionedConstant
(
    const char* const group,
    const char* entryName,
    const char* codeName,
    const unitConversion& units,
    const scalar value
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
