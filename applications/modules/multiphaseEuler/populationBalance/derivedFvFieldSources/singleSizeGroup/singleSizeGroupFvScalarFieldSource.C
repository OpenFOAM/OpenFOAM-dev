/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "singleSizeGroupFvScalarFieldSource.H"
#include "fvSource.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleSizeGroupFvScalarFieldSource::
singleSizeGroupFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    index_(dict.lookup<label>("index"))
{}


Foam::singleSizeGroupFvScalarFieldSource::
singleSizeGroupFvScalarFieldSource
(
    const singleSizeGroupFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    index_(field.index_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singleSizeGroupFvScalarFieldSource::
~singleSizeGroupFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::singleSizeGroupFvScalarFieldSource::sourceValue
(
    const fvSource& source
) const
{
    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::sizeGroup>(internalField());

    // Check the index
    const label firstIndex = fi.group().sizeGroups().first().i();
    const label lastIndex = fi.group().sizeGroups().last().i();
    if (index_ < firstIndex || index_ > lastIndex)
    {
        FatalErrorInFunction
            << "Size-group index " << index_ << " is out of range of the "
            << "indices associated with phase " << internalField().group()
            << " (" << firstIndex << " -> " << lastIndex << ")"
            << exit(FatalError);
    }

    return tmp<scalarField>(new scalarField(source.nCells(), fi.i() == index_));
}


Foam::tmp<Foam::scalarField>
Foam::singleSizeGroupFvScalarFieldSource::internalCoeff
(
    const fvSource& source
) const
{
    return tmp<scalarField>(new scalarField(source.nCells(), scalar(0)));
}


void Foam::singleSizeGroupFvScalarFieldSource::write(Ostream& os) const
{
    fvScalarFieldSource::write(os);
    writeEntry(os, "index", index_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        singleSizeGroupFvScalarFieldSource
    );
}

// ************************************************************************* //
