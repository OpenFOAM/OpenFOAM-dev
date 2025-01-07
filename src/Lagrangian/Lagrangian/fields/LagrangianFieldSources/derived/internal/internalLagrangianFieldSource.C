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

#include "internalLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::internalLagrangianFieldSource<Type>::~internalLagrangianFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::internalLagrangianFieldSource<Type>::sourceValue
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return
        LagrangianSubField<Type>::New
        (
            this->internalField().name() + ":" + source.name() + "Value",
            subMesh.sub(this->internalField())()
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::internalLagrangianFieldSource<Type>::internalCoeff
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
            dimensionedScalar(dimless, scalar(1))
        );
}


// ************************************************************************* //
