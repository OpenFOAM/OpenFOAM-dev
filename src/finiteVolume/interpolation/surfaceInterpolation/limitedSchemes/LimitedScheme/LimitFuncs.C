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

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Type>::operator()
(
    const VolField<Type>& phi
) const
{
    return Foam::magSqr(phi);
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Foam::scalar>::operator()
(
    const volScalarField& phi
) const
{
    return phi;
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Foam::symmTensor>::operator()
(
    const volSymmTensorField& phi
) const
{
    return Foam::tr(phi);
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Foam::tensor>::operator()
(
    const volTensorField& phi
) const
{
    return Foam::tr(phi);
}


template<class Type>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::rhoMagSqr<Type>::operator()
(
    const VolField<Type>& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::template lookupObject<volScalarField>("rho");
    return Foam::magSqr(phi/rho);
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::rhoMagSqr<Foam::scalar>::operator()
(
    const volScalarField& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::lookupObject<volScalarField>("rho");
    return phi/rho;
}


// ************************************************************************* //
