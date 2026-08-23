/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "surfaceThermo.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<template<class> class PrimitiveField>
void Foam::surfaceThermo::resetDerived
(
    tmp<LagrangianSubField<scalar, PrimitiveField>>& tf
) const
{
    if (tf.valid())
    {
        tf = tmp<LagrangianSubField<scalar, PrimitiveField>>
        (
            NullObjectRef<LagrangianSubField<scalar, PrimitiveField>>()
        );
    }
}


template<template<class> class PrimitiveField>
const Foam::LagrangianSubField<Foam::scalar, PrimitiveField>&
Foam::surfaceThermo::getDerived
(
    const LagrangianSubMesh& subMesh,
    tmp<LagrangianSubField<scalar, PrimitiveField>>& tf
) const
{
    if (!tf.valid())
    {
        tf =
            tmp<LagrangianSubField<scalar, PrimitiveField>>
            (
                NullObjectRef<LagrangianSubField<scalar, PrimitiveField>>()
            );
    }

    if (isNull(tf()))
    {
        correct(subMesh);
    }

    return tf();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class SurfaceThermo>
Foam::autoPtr<SurfaceThermo> Foam::surfaceThermo::New
(
    const basicLagrangianSubThermo& farThermo
)
{
    const word thermoTypeName =
        SurfaceThermo::derivedThermoName()
      + "<"
      + farThermo.mixtureName()
      + ">";

    typename SurfaceThermo::cloudConstructorTable::iterator cstrIter =
        SurfaceThermo::cloudConstructorTablePtr_->find(thermoTypeName);

    if (cstrIter == SurfaceThermo::cloudConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << SurfaceThermo::typeName << " type " << nl << nl
            << thermoTypeName << nl << nl
            << "Valid " << SurfaceThermo::typeName << " types are:" << nl
            << SurfaceThermo::cloudConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<SurfaceThermo>(cstrIter()(farThermo));
}


// ************************************************************************* //
