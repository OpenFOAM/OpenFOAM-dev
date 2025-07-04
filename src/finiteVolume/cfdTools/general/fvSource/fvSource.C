/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2025 OpenFOAM Foundation
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

#include "fvSource.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvSource, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::fvSource::infoField
(
    const word& name,
    const dimensionSet& dims,
    const scalarField& field,
    const bool print
)
{
    if (!print) return;

    Info<< indent << "min/average/max " << name
        << " = " << gMin(field) << '/' << gAverage(field) << '/' << gMax(field)
        << ' ' << dims << endl;
}


void Foam::fvSource::infoField
(
    const word& name,
    const DimensionedField<scalar, volMesh>& field,
    const bool print
)
{
    infoField(name, field.dimensions(), field.primitiveField(), print);
}


void Foam::fvSource::infoField
(
    const DimensionedField<scalar, volMesh>& field,
    const bool print
)
{
    infoField(field.name(), field, print);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvSource::fvSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvSource::~fvSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fvSource::addSupFields() const
{
    return wordList::null();
}


// ************************************************************************* //
