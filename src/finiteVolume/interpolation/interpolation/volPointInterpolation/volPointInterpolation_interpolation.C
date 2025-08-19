/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "volPointInterpolation_interpolation.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolations::volPointInterpolation<Type>::volPointInterpolation
(
    const VolField<Type>& psi
)
:
    psip_
    (
        Foam::volPointInterpolation::New(psi.mesh()).interpolate
        (
            psi,
            "volPointInterpolate(" + psi.name() + ')',
            true        // use cache
        )
    )
{}


template<class Type>
Foam::interpolations::volPointInterpolation<Type>::volPointInterpolation
(
    const volPointInterpolation<Type>& i
)
:
    psip_(i.psip_.clone())
{}


template<class Type>
Foam::interpolations::volPointInterpolation<Type>::volPointInterpolation
(
    const VolField<Type>& psi,
    tmp<PointField<Type>> psip
)
:
    psip_(psip)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::interpolations::volPointInterpolation<Type>::~volPointInterpolation()
{}



// ************************************************************************* //
