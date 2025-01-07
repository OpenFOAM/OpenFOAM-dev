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

#include "CarrierEqn.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CarrierEqn<Type>::CarrierEqn(const VolField<Type>& vf)
:
    psi(vf),
    Su(psi),
    Sp(psi)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::CarrierEqn<Type>::clear()
{
    Su.clear();
    Sp.clear();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::CarrierEqn<Type>::operator+=(const LagrangianEqn<Type>& lEqn)
{
    Su += lEqn.allDiagonalSu();
    Sp += lEqn.allDiagonalSp();
}


template<class Type>
void Foam::CarrierEqn<Type>::operator-=(const LagrangianEqn<Type>& lEqn)
{
    Su -= lEqn.allDiagonalSu();
    Sp -= lEqn.allDiagonalSp();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::operator+=
(
    fvMatrix<Type>& fvEqn,
    const CarrierEqn<Type>& cEqn
)
{
    const dimensionedScalar& deltaT = fvEqn.psi().mesh().time().deltaT();

    if (cEqn.Su.valid())
    {
        fvEqn.dimensions() -= cEqn.Su.S().dimensions()/dimTime;
        fvEqn.source() -= cEqn.Su.S().primitiveField()/deltaT.value();
    }
    if (cEqn.Sp.valid())
    {
        fvEqn.dimensions() +=
            cEqn.Sp.S().dimensions()/dimTime/fvEqn.psi().dimensions();
        fvEqn.diag() += cEqn.Sp.S().primitiveField()/deltaT.value();
    }
}


template<class Type>
void Foam::operator-=
(
    fvMatrix<Type>& fvEqn,
    const CarrierEqn<Type>& cEqn
)
{
    const dimensionedScalar& deltaT = fvEqn.psi().mesh().time().deltaT();

    if (cEqn.Su.valid())
    {
        fvEqn.dimensions() += cEqn.Su.S().dimensions()/dimTime;
        fvEqn.source() += cEqn.Su.S().primitiveField()/deltaT.value();
    }
    if (cEqn.Sp.valid())
    {
        fvEqn.dimensions() -=
            cEqn.Sp.S().dimensions()/dimTime/fvEqn.psi().dimensions();
        fvEqn.diag() -= cEqn.Sp.S().primitiveField()/deltaT.value();
    }
}


// ************************************************************************* //
