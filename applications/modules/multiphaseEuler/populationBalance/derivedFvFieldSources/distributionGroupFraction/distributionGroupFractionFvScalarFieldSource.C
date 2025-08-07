/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "distributionGroupFractionFvScalarFieldSource.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::distributionGroupFractionFvScalarFieldSource::eta
(
    const fvSource& model
) const
{
    if (!etaPtr_.valid())
    {
        etaPtr_.set
        (
            new scalar
            (
                popBal().etaV(i(), distribution_()).value()
            )
        );

        if (debug)
        {
            Info<< typeName << ": Group #" << i() << " is receiving "
                << 100*etaPtr_() << "\% of source " << model.name() << endl;
        }
    }

    return etaPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionGroupFractionFvScalarFieldSource::
distributionGroupFractionFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    groupPropertyFvScalarField(iF),
    distribution_
    (
        distribution::New(dimLength, dict.subDict("distribution"), 3, -1)
    ),
    etaPtr_(nullptr)
{}


Foam::distributionGroupFractionFvScalarFieldSource::
distributionGroupFractionFvScalarFieldSource
(
    const distributionGroupFractionFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    groupPropertyFvScalarField(iF),
    distribution_(field.distribution_, false),
    etaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionGroupFractionFvScalarFieldSource::
~distributionGroupFractionFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::distributionGroupFractionFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        DimensionedField<scalar, volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "SourceValue",
            this->internalField().mesh(),
            dimensionedScalar(dimless, eta(model))
        );
}


Foam::tmp<Foam::scalarField>
Foam::distributionGroupFractionFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return tmp<scalarField>(new scalarField(source.size(), eta(model)));
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::distributionGroupFractionFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        DimensionedField<scalar, volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "InternalCoeff",
            this->internalField().mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


Foam::tmp<Foam::scalarField>
Foam::distributionGroupFractionFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return tmp<scalarField>(new scalarField(source.size(), scalar(0)));
}


void Foam::distributionGroupFractionFvScalarFieldSource::write
(
    Ostream& os
) const
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
        distributionGroupFractionFvScalarFieldSource
    );

    // Backwards compatible lookup as "distributionSizeGroup"
    addNamedToRunTimeSelectionTable
    (
        fvScalarFieldSource,
        distributionGroupFractionFvScalarFieldSource,
        dictionary,
        distributionSizeGroup
    );
}

// ************************************************************************* //
