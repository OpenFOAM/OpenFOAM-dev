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

#include "nucleationInterfacialCurvatureFvScalarFieldSource.H"
#include "fvSource.H"
#include "nucleation.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nucleationInterfacialCurvatureFvScalarFieldSource::
nucleationInterfacialCurvatureFvScalarFieldSource
(
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fvScalarFieldSource(iF, dict)
{}


Foam::nucleationInterfacialCurvatureFvScalarFieldSource::
nucleationInterfacialCurvatureFvScalarFieldSource
(
    const nucleationInterfacialCurvatureFvScalarFieldSource& field,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvScalarFieldSource(field, iF)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nucleationInterfacialCurvatureFvScalarFieldSource::
~nucleationInterfacialCurvatureFvScalarFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationInterfacialCurvatureFvScalarFieldSource::sourceValue
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    NotImplemented;
    return tmp<DimensionedField<scalar, volMesh>>(nullptr);
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationInterfacialCurvatureFvScalarFieldSource::internalCoeff
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
            dimensionedScalar(dimless, scalar(-1))
        );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::nucleationInterfacialCurvatureFvScalarFieldSource::sourceCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return 6/refCast<const fv::nucleation>(model).d();
}


void Foam::nucleationInterfacialCurvatureFvScalarFieldSource::write
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
        nucleationInterfacialCurvatureFvScalarFieldSource
    );
}

// ************************************************************************* //
