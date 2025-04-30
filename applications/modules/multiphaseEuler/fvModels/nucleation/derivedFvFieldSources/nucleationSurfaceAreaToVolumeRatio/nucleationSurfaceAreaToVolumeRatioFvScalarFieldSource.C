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

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    // Determine the size-group index from the end of the field name
    const word member = this->internalField().member();
    string::size_type memberChari = member.size();
    while (memberChari && isdigit(member[memberChari - 1])) memberChari --;
    const label i = atoi(member(memberChari, member.size()).c_str());

    // Access the size-group
    const diameterModels::sizeGroup& fi =
        refCast<const diameterModels::velocityGroup>
        (
            db()
           .lookupObject<phaseSystem>(phaseSystem::propertiesName)
           .phases()[this->internalField().group()]
           .diameter()
        ).popBal().sizeGroups()[i];

    tmp<volScalarField::Internal> d =
        refCast<const fv::nucleation>(model).d();

    return
        fi.group().popBal().etaV
        (
            fi.i(),
            constant::mathematical::pi/6*pow3(d)
        );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationSurfaceAreaToVolumeRatioFvScalarFieldSource::internalCoeff
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
