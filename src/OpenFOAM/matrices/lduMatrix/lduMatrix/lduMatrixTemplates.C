/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Description
    lduMatrix member H operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::lduMatrix::H(const Field<Type>& psi) const
{
    tmp<Field<Type>> tHpsi
    (
        new Field<Type>(lduAddr().size(), Zero)
    );

    if (lowerPtr_ || upperPtr_)
    {
        Field<Type> & Hpsi = tHpsi.ref();

        Type* __restrict__ HpsiPtr = Hpsi.begin();

        const Type* __restrict__ psiPtr = psi.begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* __restrict__ lowerPtr = lower().begin();
        const scalar* __restrict__ upperPtr = upper().begin();

        const label nFaces = upper().size();

        for (label face=0; face<nFaces; face++)
        {
            HpsiPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
            HpsiPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
        }
    }

    return tHpsi;
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::lduMatrix::H(const tmp<Field<Type>>& tpsi) const
{
    tmp<Field<Type>> tHpsi(H(tpsi()));
    tpsi.clear();
    return tHpsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::lduMatrix::faceH(const Field<Type>& psi) const
{
    if (lowerPtr_ || upperPtr_)
    {
        const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
        const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

        const labelUList& l = lduAddr().lowerAddr();
        const labelUList& u = lduAddr().upperAddr();

        tmp<Field<Type>> tfaceHpsi(new Field<Type> (Lower.size()));
        Field<Type> & faceHpsi = tfaceHpsi.ref();

        for (label face=0; face<l.size(); face++)
        {
            faceHpsi[face] =
                Upper[face]*psi[u[face]]
              - Lower[face]*psi[l[face]];
        }

        return tfaceHpsi;
    }
    else
    {
        FatalErrorInFunction
            << "Cannot calculate faceH"
               " the matrix does not have any off-diagonal coefficients."
            << exit(FatalError);

        return tmp<Field<Type>>(nullptr);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::lduMatrix::faceH(const tmp<Field<Type>>& tpsi) const
{
    tmp<Field<Type>> tfaceHpsi(faceH(tpsi()));
    tpsi.clear();
    return tfaceHpsi;
}


// ************************************************************************* //
