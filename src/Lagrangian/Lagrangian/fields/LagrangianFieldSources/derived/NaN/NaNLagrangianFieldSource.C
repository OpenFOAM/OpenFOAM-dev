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

#include "NaNLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::NaNLagrangianFieldSource<Type>::
~NaNLagrangianFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::NaNLagrangianFieldSource<Type>::sourceValue
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return
        LagrangianSubField<Type>::New
        (
            this->internalField().name() + ":" + source.name() + "Coeff",
            subMesh,
            dimensioned<Type>
            (
                "NaN",
                this->internalDimensions(),
                pTraits<Type>::nan
            )
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::NaNLagrangianFieldSource<Type>::internalCoeff
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
            dimensionedScalar
            (
                "NaN",
                dimless,
                pTraits<scalar>::nan
            )
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::NaNLagrangianFieldSource<Type>::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return
        LagrangianSubField<Type>::New
        (
            this->internalField().name() + ":" + injection.name() + "Coeff",
            subMesh,
            dimensioned<Type>
            (
                "NaN",
                this->internalDimensions(),
                pTraits<Type>::nan
            )
        );
}


// ************************************************************************* //
