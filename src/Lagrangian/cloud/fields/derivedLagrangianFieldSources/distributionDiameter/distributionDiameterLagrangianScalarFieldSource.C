/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "distributionDiameterLagrangianScalarFieldSource.H"
#include "uniformSizeNumberLagrangianScalarFieldSource.H"
#include "LagrangianFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionDiameterLagrangianScalarFieldSource::
distributionDiameterLagrangianScalarFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianScalarFieldSource(iIo, dict),
    distribution_
    (
        distribution::New
        (
            internalDimensions(),
            dict.subDict("distribution"),
            0,
            randomGenerator::seed(iIo.name() + ':' + dict.dictName())
        )
    ),
    timeIndex_(-1),
    numberName_(dict.lookupOrDefault<word>("number", "number"))
{}


Foam::distributionDiameterLagrangianScalarFieldSource::
distributionDiameterLagrangianScalarFieldSource
(
    const distributionDiameterLagrangianScalarFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianScalarFieldSource(field, iIo),
    distribution_(field.distribution_, false),
    timeIndex_(-1),
    numberName_(field.numberName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionDiameterLagrangianScalarFieldSource::
~distributionDiameterLagrangianScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::distributionDiameterLagrangianScalarFieldSource::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    // Look up the distribution number source
    const uniformSizeNumberLagrangianScalarFieldSource& numberUsnFs =
        fieldSourceCast<scalar, uniformSizeNumberLagrangianScalarFieldSource>
        (
            numberName_,
            injection
        );

    // The sample size exponent of the distribution is obtained from the number
    // source. This wasn't available during construction, so it is set here if
    // this is the first execution.
    if (timeIndex_ == -1)
    {
        distribution_ = distribution::New(distribution_, numberUsnFs.sampleQ());
    }

    // Restart the distribution if the time index has not changed. This ensures
    // multiple evaluations from different conditions and/or iterations
    // generate the same values
    distribution_->start(timeIndex_ == db().time().timeIndex());
    timeIndex_ = db().time().timeIndex();

    // Sample the distribution and return the result as a sub-field
    return
        LagrangianSubScalarField::New
        (
            internalField().name() + ":" + injection.name(),
            subMesh,
            internalDimensions(),
            distribution_->sample(subMesh.size())
        );
}


void Foam::distributionDiameterLagrangianScalarFieldSource::write
(
    Ostream& os
) const
{
    LagrangianScalarFieldSource::write(os);

    writeEntry(os, "distribution", internalDimensions(), distribution_());
    writeEntryIfDifferent<word>(os, "number", "number", numberName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLagrangianTypeFieldSource
    (
        LagrangianScalarFieldSource,
        distributionDiameterLagrangianScalarFieldSource
    );
}

// ************************************************************************* //
