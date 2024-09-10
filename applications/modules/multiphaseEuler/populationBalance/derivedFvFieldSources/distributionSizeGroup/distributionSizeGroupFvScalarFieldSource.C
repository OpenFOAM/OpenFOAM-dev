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

#include "distributionSizeGroupFvScalarFieldSource.H"
#include "fvSource.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionSizeGroupFvScalarFieldSource::
distributionSizeGroupFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    distribution_
    (
        distribution::New(dimLength, dict.subDict("distribution"), 3, -1)
    ),
    etaPtr_(nullptr)
{}


Foam::distributionSizeGroupFvScalarFieldSource::
distributionSizeGroupFvScalarFieldSource
(
    const distributionSizeGroupFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    distribution_(field.distribution_, false),
    etaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionSizeGroupFvScalarFieldSource::
~distributionSizeGroupFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::distributionSizeGroupFvScalarFieldSource::sourceValue
(
    const fvSource& source
) const
{
    if (!etaPtr_.valid())
    {
        const diameterModels::sizeGroup& fi =
            refCast<const diameterModels::sizeGroup>(internalField());

        etaPtr_.set
        (
            new scalar
            (
                fi.group().popBal().etaV(fi.i(), distribution_()).value()
            )
        );

        if (debug)
        {
            Info<< typeName << ": Size group #" << fi.i() << " is receiving "
                << 100*etaPtr_() << "\% of source " << source.name() << endl;
        }
    }

    return tmp<scalarField>(new scalarField(source.nCells(), etaPtr_()));
}


Foam::tmp<Foam::scalarField>
Foam::distributionSizeGroupFvScalarFieldSource::internalCoeff
(
    const fvSource& source
) const
{
    return tmp<scalarField>(new scalarField(source.nCells(), scalar(0)));
}


void Foam::distributionSizeGroupFvScalarFieldSource::write(Ostream& os) const
{
    fvScalarFieldSource::write(os);
    writeEntry(os, "distribution", dimLength, distribution_(), true, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        distributionSizeGroupFvScalarFieldSource
    );
}

// ************************************************************************* //
