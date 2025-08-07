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

#include "nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource.H"
#include "populationBalanceModel.H"
#include "nucleation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource::
nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    groupPropertyFvScalarField(iF)
{}


Foam::nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource::
nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource
(
    const nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    groupPropertyFvScalarField(iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource::
~nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    // Get the surface area volume ratio of the nucleates
    tmp<volScalarField::Internal> kappa =
        6/refCast<const fv::nucleation>(model).d();

    // Scale the value by the source for the corresponding group fraction
    return popBal().f(i()).sources()[model.name()].value(model, source)*kappa;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource::internalCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    // Nucleation is always an "inflow" to the nucleating phase, so the source
    // should be fully explicit
    return
        DimensionedField<scalar, volMesh>::New
        (
            model.name() + ":" + this->internalField().name() + "InternalCoeff",
            this->internalField().mesh(),
            dimensionedScalar(dimless, scalar(0))
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource
    );

    // Backwards compatible lookup as "nucleationSurfaceAreaVolumeRatio"
    addNamedToRunTimeSelectionTable
    (
        fvScalarFieldSource,
        nucleationGroupSurfaceAreaVolumeRatioFvScalarFieldSource,
        dictionary,
        nucleationSurfaceAreaVolumeRatio
    );
}

// ************************************************************************* //
