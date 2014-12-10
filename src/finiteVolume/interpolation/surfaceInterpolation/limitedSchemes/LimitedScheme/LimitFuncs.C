/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

namespace Foam
{

namespace limitFuncs
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline tmp<volScalarField> magSqr<Type>::operator()
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    return Foam::magSqr(phi);
}


template<>
inline tmp<volScalarField> magSqr<scalar>::operator()
(
    const volScalarField& phi
) const
{
    return phi;
}


template<>
inline tmp<volScalarField> magSqr<tensor>::operator()
(
    const volTensorField& phi
) const
{
    return Foam::tr(phi);
}


template<class Type>
inline tmp<volScalarField> rhoMagSqr<Type>::operator()
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::template lookupObject<volScalarField>("rho");
    return Foam::magSqr(phi/rho);
}


template<>
inline tmp<volScalarField> rhoMagSqr<scalar>::operator()
(
    const volScalarField& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::lookupObject<volScalarField>("rho");
    return phi/rho;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace limitFuncs
} // End namespace Foam

// ************************************************************************* //
