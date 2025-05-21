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

#include "nucleationSurfaceAreaVolumeRatioFvScalarFieldSource.H"
#include "volFields.H"
#include "fvSource.H"
#include "sizeGroup.H"
#include "nucleation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationSurfaceAreaVolumeRatioFvScalarFieldSource::
nucleationSurfaceAreaVolumeRatioFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    secondaryPropertyFvScalarFieldSource(iF)
{}


Foam::nucleationSurfaceAreaVolumeRatioFvScalarFieldSource::
nucleationSurfaceAreaVolumeRatioFvScalarFieldSource
(
    const nucleationSurfaceAreaVolumeRatioFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    secondaryPropertyFvScalarFieldSource(iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationSurfaceAreaVolumeRatioFvScalarFieldSource::
~nucleationSurfaceAreaVolumeRatioFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationSurfaceAreaVolumeRatioFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    // Get the diameter of the nucleates
    tmp<volScalarField::Internal> d = refCast<const fv::nucleation>(model).d();

    // Return the ratio, scaled by the source for the corresponding size-group
    return fi().sources()[model.name()].value(model, source)*6/d;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationSurfaceAreaVolumeRatioFvScalarFieldSource::internalCoeff
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
        nucleationSurfaceAreaVolumeRatioFvScalarFieldSource
    );
}

// ************************************************************************* //
