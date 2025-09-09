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
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class PsiFieldType>
Foam::tmp<Foam::VolInternalField<Type>>
Foam::CarrierEqn<Type>::residual
(
    const dimensionSet& dims,
    const word& psiName,
    const PsiFieldType& psi
) const
{
    // Check the dimensions
    if (Su.valid())
    {
        Su.S().dimensions() = dims*dimTime*dimVolume;
    }
    if (Sp.valid())
    {
        Sp.S().dimensions()*psi.dimensions() = dims*dimTime*dimVolume;
    }

    // Build the residual
    if (Su.valid() && Sp.valid())
    {
        return (Su.S() + Sp.S()*psi)/mesh_.time().deltaT()/mesh_.V();
    }
    else if (Su.valid())
    {
        return Su.S()/mesh_.time().deltaT()/mesh_.V();
    }
    else if (Sp.valid())
    {
        return Sp.S()*psi/mesh_.time().deltaT()/mesh_.V();
    }
    else
    {
        return
            VolInternalField<Type>::New
            (
                name_ + "&" + psiName,
                mesh_,
                dimensionedScalar(dims, scalar(0))
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CarrierEqn<Type>::CarrierEqn(const VolField<Type>& psi)
:
    CarrierEqnBase
    (
        IOobject::groupName
        (
            IOobject::member(psi.name()) + "Eqn",
            IOobject::group(psi.name())
        ),
        psi.mesh()
    ),
    psi_(psi),
    Su(*this),
    Sp(*this)
{}


template<class Type>
Foam::CarrierEqn<Type>::CarrierEqn
(
    const word& psiName,
    const fvMesh& mesh
)
:
    CarrierEqnBase
    (
        IOobject::groupName
        (
            IOobject::member(psiName) + "Eqn",
            IOobject::group(psiName)
        ),
        mesh
    ),
    psi_(NullObjectRef<VolField<Type>>()),
    Su(*this),
    Sp(*this)
{}


template<class Type>
Foam::CarrierEqn<Type>::CarrierEqn(const CarrierEqn<Type>& eqn)
:
    CarrierEqnBase(eqn.name(), eqn.mesh()),
    psi_(eqn.psi_),
    Su(eqn.Su),
    Sp(eqn.Sp)
{}


template<class Type>
Foam::CarrierEqn<Type>::CarrierEqn(const tmp<CarrierEqn<Type>>& tEqn)
:
    CarrierEqnBase(tEqn().name(), tEqn().mesh()),
    psi_(tEqn().psi_),
    Su
    (
        tEqn.isTmp()
      ? CarrierCoeff<Type, false>(tEqn.ref().Su, true)
      : CarrierCoeff<Type, false>(tEqn().Su)
    ),
    Sp
    (
        tEqn.isTmp()
      ? CarrierCoeff<scalar, true>(tEqn.ref().Sp, true)
      : CarrierCoeff<scalar, true>(tEqn().Sp)
    )
{
    tEqn.clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::CarrierEqn<Type>::clear()
{
    Su.clear();
    Sp.clear();
}


template<class Type>
Foam::tmp<Foam::VolInternalField<Type>>
Foam::CarrierEqn<Type>::residual
(
    const dimensionSet& dims
) const
{
    return
        notNull(psi_)
      ? residual(dims, psi_.name(), psi_)
      : residual(dims, "1", geometricOneField());
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::CarrierEqn<Type>::operator+=(const LagrangianEqn<Type>& lEqn)
{
    Su += lEqn.allDiagonalSu();
    Sp += lEqn.allDiagonalSp();
}


template<class Type>
void Foam::CarrierEqn<Type>::operator+=(const tmp<LagrangianEqn<Type>>& tlEqn)
{
    operator+=(tlEqn());
    tlEqn.clear();
}


template<class Type>
void Foam::CarrierEqn<Type>::operator-=(const LagrangianEqn<Type>& lEqn)
{
    Su -= lEqn.allDiagonalSu();
    Sp -= lEqn.allDiagonalSp();
}


template<class Type>
void Foam::CarrierEqn<Type>::operator-=(const tmp<LagrangianEqn<Type>>& tlEqn)
{
    operator-=(tlEqn());
    tlEqn.clear();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

#define CARRIER_EQN_SCALAR_OPERATOR(Op, EqOp)                                  \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::CarrierEqn<Type>> Foam::operator Op                        \
    (                                                                          \
        const CarrierEqn<Type>& eqn,                                           \
        const dimensionedScalar& s                                             \
    )                                                                          \
    {                                                                          \
        return tmp<CarrierEqn<Type>>(eqn) Op s;                                \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::CarrierEqn<Type>> Foam::operator Op                        \
    (                                                                          \
        const tmp<CarrierEqn<Type>>& tEqn,                                     \
        const dimensionedScalar& s                                             \
    )                                                                          \
    {                                                                          \
        tmp<CarrierEqn<Type>> tResult(new CarrierEqn<Type>(tEqn));             \
                                                                               \
        tResult.ref().Su EqOp s;                                               \
        tResult.ref().Sp EqOp s;                                               \
                                                                               \
        return tResult;                                                        \
    }

#define CARRIER_COMMUTATIVE_SCALAR_EQN_OPERATOR(Op, EqOp)                      \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::CarrierEqn<Type>> Foam::operator Op                        \
    (                                                                          \
        const dimensionedScalar& s,                                            \
        const CarrierEqn<Type>& eqn                                            \
    )                                                                          \
    {                                                                          \
        return eqn Op s;                                                       \
    }                                                                          \
                                                                               \
    template<class Type>                                                       \
    Foam::tmp<Foam::CarrierEqn<Type>> Foam::operator Op                        \
    (                                                                          \
        const dimensionedScalar& s,                                            \
        const tmp<CarrierEqn<Type>>& tEqn                                      \
    )                                                                          \
    {                                                                          \
        return tEqn Op s;                                                      \
    }

//- Multiplication
CARRIER_EQN_SCALAR_OPERATOR(*, *=)
CARRIER_COMMUTATIVE_SCALAR_EQN_OPERATOR(*, *=)

//- Division
CARRIER_EQN_SCALAR_OPERATOR(/, /=)

#undef CARRIER_EQN_SCALAR_OPERATOR
#undef CARRIER_COMMUTATIVE_SCALAR_EQN_OPERATOR

template<class Type>
void Foam::operator+=(fvMatrix<Type>& fvEqn, const CarrierEqn<Type>& cEqn)
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
            cEqn.Sp.S().dimensions()/dimTime*fvEqn.psi().dimensions();
        fvEqn.diag() += cEqn.Sp.S().primitiveField()/deltaT.value();
    }
}


template<class Type>
void Foam::operator+=(fvMatrix<Type>& fvEqn, const tmp<CarrierEqn<Type>>& tcEqn)
{
    fvEqn += tcEqn();
    tcEqn.clear();
}


template<class Type>
void Foam::operator-=(fvMatrix<Type>& fvEqn, const CarrierEqn<Type>& cEqn)
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
            cEqn.Sp.S().dimensions()/dimTime*fvEqn.psi().dimensions();
        fvEqn.diag() -= cEqn.Sp.S().primitiveField()/deltaT.value();
    }
}


template<class Type>
void Foam::operator-=(fvMatrix<Type>& fvEqn, const tmp<CarrierEqn<Type>>& tcEqn)
{
    fvEqn -= tcEqn();
    tcEqn.clear();
}


// ************************************************************************* //
