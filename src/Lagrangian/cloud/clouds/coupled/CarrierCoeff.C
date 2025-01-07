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

#include "CarrierCoeff.H"
#include "LagrangiancAccumulate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, bool Implicit>
void Foam::CarrierCoeff<Type, Implicit>::initialise
(
    const LagrangianSubField<Type>& lField
)
{
    if (valid()) return;

    S_.set
    (
        new DimensionedField<Type, volMesh>
        (
            IOobject
            (
                lField.mesh().mesh().name() + ":" + lField.name(),
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensioned<Type>(lField.dimensions(), Zero)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, bool Implicit>
template<class PsiType>
Foam::CarrierCoeff<Type, Implicit>::CarrierCoeff(const VolField<PsiType>& psi)
:
    mesh_(psi.mesh())
{}


template<class Type, bool Implicit>
Foam::CarrierCoeff<Type, Implicit>::CarrierCoeff
(
    const CarrierCoeff<Type, Implicit>& coeff
)
:
    mesh_(coeff.mesh_),
    S_(coeff.S_, false)
{}


template<class Type, bool Implicit>
Foam::CarrierCoeff<Type, Implicit>::CarrierCoeff
(
    CarrierCoeff<Type, Implicit>& coeff,
    const bool reuse
)
:
    mesh_(coeff.mesh_),
    S_(coeff.S_, reuse)
{}


template<class Type, bool Implicit>
Foam::CarrierCoeff<Type, Implicit>::CarrierCoeff
(
    CarrierCoeff<Type, Implicit>&& coeff
)
:
    mesh_(coeff.mesh_),
    S_(coeff.S_, true)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, bool Implicit>
bool Foam::CarrierCoeff<Type, Implicit>::valid() const
{
    return S_.valid();
}


template<class Type, bool Implicit>
void Foam::CarrierCoeff<Type, Implicit>::clear()
{
    S_.clear();
}


template<class Type, bool Implicit>
const Foam::DimensionedField<Type, Foam::volMesh>&
Foam::CarrierCoeff<Type, Implicit>::S() const
{
    return S_();
}


template<class Type, bool Implicit>
void Foam::CarrierCoeff<Type, Implicit>::negate()
{
    if (valid()) S_() = -S_();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type, bool Implicit>
void Foam::CarrierCoeff<Type, Implicit>::operator+=
(
    const LagrangianCoeff<Type, Implicit>& lCoeff
)
{
    if (!lCoeff.valid()) return;

    initialise(lCoeff.S());

    Lagrangianc::accumulate(lCoeff.S(), S_(), lCoeff.eqn().name());
}


template<class Type, bool Implicit>
void Foam::CarrierCoeff<Type, Implicit>::operator-=
(
    const LagrangianCoeff<Type, Implicit>& lCoeff
)
{
    if (!lCoeff.valid()) return;

    initialise(lCoeff.S());

    Lagrangianc::accumulate((-lCoeff.S())(), S_(), lCoeff.eqn().name());
}


// ************************************************************************* //
