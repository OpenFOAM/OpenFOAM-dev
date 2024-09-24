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

#include "nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource.H"
#include "fvSource.H"
#include "populationBalanceModel.H"
#include "nucleation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::
nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict)
{}


Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::
nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource
(
    const nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::
~nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::sourceValue
(
    const fvSource& source
) const
{
    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::sizeGroup>(this->internalField());

    tmp<volScalarField::Internal> d =
        refCast<const fv::nucleation>(source).d();

    tmp<volScalarField::Internal> eta =
        fi.group().popBal().etaV(fi.i(), constant::mathematical::pi/6*pow3(d));

    return tmp<scalarField>(new scalarField(eta*6/d, source.cells()));
}


Foam::tmp<Foam::scalarField>
Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::internalCoeff
(
    const fvSource& source
) const
{
    // Nucleation is always an "inflow" to the nucleating phase, so the source
    // should be fully explicit
    return tmp<scalarField>(new scalarField(source.nCells(), scalar(0)));
}


void Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::write
(
    Ostream& os
) const
{
    fvScalarFieldSource::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource
    );
}

// ************************************************************************* //
