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

#include "LagrangianSp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianSp<Type>::LagrangianSp(const LagrangianSp<Type>& Sp)
:
    LagrangianCoeff<scalar, true>(Sp)
{}


template<class Type>
Foam::LagrangianSp<Type>::LagrangianSp(LagrangianSp<Type>& Sp, const bool reuse)
:
    LagrangianCoeff<scalar, true>(Sp, reuse)
{}


template<class Type>
Foam::LagrangianSp<Type>::LagrangianSp(LagrangianSp<Type>&& Sp)
:
    LagrangianCoeff<scalar, true>(move(Sp))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Foam::scalar, true>>
Foam::LagrangianSp<Type>::A() const
{
    return *this;
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianSp<Type>::H() const
{
    return
        tmp<LagrangianCoeff<Type, false>>
        (
            new LagrangianCoeff<Type, false>(eqn())
        );
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianSp<Type>::Su() const
{
    return Su(static_cast<const LagrangianEqn<Type>&>(this->eqn()));
}


template<class Type>
Foam::tmp<Foam::LagrangianCoeff<Type, false>>
Foam::LagrangianSp<Type>::Su(const LagrangianEqn<Type>& eqn) const
{
    return
        tmp<LagrangianCoeff<Type, false>>
        (
            valid()
          ? new LagrangianCoeff<Type, false>(eqn, S()*eqn.psi())
          : new LagrangianCoeff<Type, false>(eqn)
        );
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianSp<Type>::operator+=(const LagrangianSp<Type>& Sp)
{
    static_cast<LagrangianCoeff<scalar, true>&>(*this) += Sp;
}


template<class Type>
void Foam::LagrangianSp<Type>::operator-=(const LagrangianSp<Type>& Sp)
{
    static_cast<LagrangianCoeff<scalar, true>&>(*this) -= Sp;
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>> Foam::operator/
(
    const LagrangianCoeff<Type, false>& Su,
    const LagrangianSp<Type>& Sp
)
{
    return Su.S()/Sp.S();
}


// ************************************************************************* //
