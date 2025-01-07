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

#include "Function1LagrangianFieldSource.H"
#include "LagrangianMesh.H"
#include "LagrangianModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
template<class Derived>
Foam::Function1LagrangianFieldSource<Type>::Function1LagrangianFieldSource
(
    const Derived& field
)
:
    field_(static_cast<const LagrangianFieldSource<Type>&>(field))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
template<class OtherType>
Foam::tmp<Foam::LagrangianSubField<OtherType>>
Foam::Function1LagrangianFieldSource<Type>::value
(
    const LagrangianModel& model,
    const LagrangianSubMesh& subMesh,
    const dimensionSet& dims,
    const Function1<OtherType>& function
)
{
    const objectRegistry& db = subMesh.mesh();

    const scalar t1 = db.time().value();
    const scalar t0 = t1 - db.time().deltaTValue();

    const word name = model.name() + ":" + function.name();

    if
    (
        db.template foundObject<LagrangianScalarInternalDynamicField>
        (
            LagrangianMesh::fractionName
        )
    )
    {
        const LagrangianSubScalarSubField fractionSf
        (
            subMesh.sub
            (
                db.template lookupObject<LagrangianScalarInternalDynamicField>
                (
                    LagrangianMesh::fractionName
                )
            )
        );

        return
            LagrangianSubField<OtherType>::New
            (
                name,
                subMesh,
                dims,
                function.value(t0 + fractionSf*(t1 - t0))
            );
    }
    else
    {
        return
            LagrangianSubField<OtherType>::New
            (
                name,
                subMesh,
                dimensioned<OtherType>(dims, function.value(t0))
            );
    }
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::Function1LagrangianFieldSource<Type>::value
(
    const LagrangianModel& model,
    const LagrangianSubMesh& subMesh,
    const Function1<Type>& function
) const
{
    return value(model, subMesh, field_.internalDimensions(), function);
}


// ************************************************************************* //
