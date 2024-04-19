/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "steadyStateD2dt2Scheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
steadyStateD2dt2Scheme<Type>::fvcD2dt2
(
    const VolField<Type>& vf
)
{
    return VolField<Type>::New
    (
        "d2dt2("+vf.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            vf.dimensions()/dimTime/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<VolField<Type>>
steadyStateD2dt2Scheme<Type>::fvcD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return VolField<Type>::New
    (
        "d2dt2("+rho.name()+','+vf.name()+')',
        mesh(),
        dimensioned<Type>
        (
            "0",
            rho.dimensions()*vf.dimensions()/dimTime/dimTime,
            Zero
        )
    );
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateD2dt2Scheme<Type>::fvmD2dt2
(
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVolume/dimTime/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateD2dt2Scheme<Type>::fvmD2dt2
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVolume/dimTime/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateD2dt2Scheme<Type>::fvmD2dt2
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVolume/dimTime/dimTime
        )
    );

    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
