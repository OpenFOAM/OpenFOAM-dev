/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "fvcMagSqrGradGrad.H"
#include "fvcGrad.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<volScalarField> magSqrGradGrad
(
    const VolField<Type>& vf
)
{
    tmp<volScalarField> tMagSqrGradGrad
    (
        magSqr(fvc::grad(fvc::grad(vf.component(0))))
    );

    // Loop over other vector field components
    for (direction cmpt = 1; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        tMagSqrGradGrad.ref() +=
            magSqr(fvc::grad(fvc::grad(vf.component(cmpt))))();
    }

    return tMagSqrGradGrad;
}


template<class Type>
tmp<volScalarField>
magSqrGradGrad
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<volScalarField> tMagSqrGradGrad(fvc::magSqrGradGrad(tvf()));
    tvf.clear();
    return tMagSqrGradGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
