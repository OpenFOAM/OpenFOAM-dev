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

#include "uniformFixedValueLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueLagrangianFieldSource<Type>::
uniformFixedValueLagrangianFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianFieldSource<Type>(iIo, dict),
    Function1LagrangianFieldSource<Type>(*this),
    uniformValue_
    (
        Function1<Type>::New
        (
            "uniformValue",
            iIo.time().userUnits(),
            this->internalDimensions(),
            dict
        )
    )
{}


template<class Type>
Foam::uniformFixedValueLagrangianFieldSource<Type>::
uniformFixedValueLagrangianFieldSource
(
    const uniformFixedValueLagrangianFieldSource<Type>& field,
    const regIOobject& iIo
)
:
    LagrangianFieldSource<Type>(field, iIo),
    Function1LagrangianFieldSource<Type>(*this),
    uniformValue_(field.uniformValue_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueLagrangianFieldSource<Type>::
~uniformFixedValueLagrangianFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::uniformFixedValueLagrangianFieldSource<Type>::sourceValue
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return value(source, subMesh, uniformValue_());
}


template<class Type>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::uniformFixedValueLagrangianFieldSource<Type>::internalCoeff
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return
        LagrangianSubScalarField::New
        (
            this->internalField().name() + ":" + source.name() + "Coeff",
            subMesh,
            dimensionedScalar(dimless, Zero)
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::uniformFixedValueLagrangianFieldSource<Type>::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return value(injection, subMesh, uniformValue_());
}


template<class Type>
void Foam::uniformFixedValueLagrangianFieldSource<Type>::write
(
    Ostream& os
) const
{
    LagrangianFieldSource<Type>::write(os);

    writeEntry(os, uniformValue_());
}


// ************************************************************************* //
