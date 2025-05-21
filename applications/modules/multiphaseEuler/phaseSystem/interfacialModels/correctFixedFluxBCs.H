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

Namespace
    Foam

Description
    Function for correcting model coefficients on patches with fixed flux
    boundary conditions

\*---------------------------------------------------------------------------*/

#ifndef correctFixedFluxBCs_H
#define correctFixedFluxBCs_H

#include "phaseInterface.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
tmp<GeometricField<Type, GeoMesh>> correctFixedFluxBCs
(
    const phaseInterface& interface,
    tmp<GeometricField<Type, GeoMesh>> tfield
)
{
    typedef GeometricField<Type, GeoMesh> typeGeoField;

    typename typeGeoField::Boundary& fieldBf = tfield.ref().boundaryFieldRef();

    forAll(fieldBf, patchi)
    {
        if
        (
            (
                !interface.phase1().stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    interface.phase1().phi()().boundaryField()[patchi]
                )
            )
         || (
                !interface.phase2().stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    interface.phase2().phi()().boundaryField()[patchi]
                )
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }

    return tfield;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
