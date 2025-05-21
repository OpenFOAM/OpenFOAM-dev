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

#include "uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::
uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict),
    secondaryPropertyFvScalarFieldSource(iF),
    uniformValue_
    (
        Function1<scalar>::New
        (
            "uniformValue",
            this->db().time().userUnits(),
            iF.dimensions(),
            dict
        )
    )
{}


Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::
uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource
(
    const uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF),
    secondaryPropertyFvScalarFieldSource(iF),
    uniformValue_(field.uniformValue_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::
~uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    // Return the ratio, scaled by the source for the corresponding size-group
    return
        fi().sources()[model.name()].value(model, source)
       *dimensionedScalar
        (
            internalField().dimensions(),
            uniformValue_->value(this->db().time().value())
        );
}


Foam::tmp<Foam::scalarField>
Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    // Return the ratio, scaled by the source for the corresponding size-group
    return
        fi().sources()[model.name()].value(model, source, cells)
       *uniformValue_->value(this->db().time().value());
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::
internalCoeff
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
Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::
internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return tmp<scalarField>(new scalarField(source.size(), scalar(0)));
}


void Foam::uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
    writeEntry
    (
        os,
        this->db().time().userUnits(),
        this->internalField().dimensions(),
        uniformValue_()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        uniformFixedValueSurfaceAreaVolumeRatioFvScalarFieldSource
    );
}

// ************************************************************************* //
