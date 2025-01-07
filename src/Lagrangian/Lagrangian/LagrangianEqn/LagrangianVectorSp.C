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
#include "LagrangianEqn.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LagrangianSp<Foam::vector>::makeTensor()
{
    if (tensorCoeff_.valid()) return;

    if (!scalarCoeff_.valid()) return;

    tensorCoeff_ += tensor::I*scalarCoeff_.S();

    scalarCoeff_ *= Zero;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianSp<Foam::vector>::LagrangianSp(const LagrangianEqn<vector>& eqn)
:
    scalarCoeff_(eqn),
    tensorCoeff_(eqn)
{}


Foam::LagrangianSp<Foam::vector>::LagrangianSp
(
    const LagrangianSp<vector>& coeff
)
:
    tmp<LagrangianSp<vector>>::refCount(),
    scalarCoeff_(coeff.scalarCoeff_),
    tensorCoeff_(coeff.tensorCoeff_)
{}


Foam::LagrangianSp<Foam::vector>::LagrangianSp
(
    LagrangianSp<vector>& coeff,
    const bool reuse
)
:
    scalarCoeff_(coeff.scalarCoeff_, reuse),
    tensorCoeff_(coeff.tensorCoeff_, reuse)
{}


Foam::LagrangianSp<Foam::vector>::LagrangianSp(LagrangianSp<vector>&& coeff)
:
    tmp<LagrangianSp<vector>>::refCount(),
    scalarCoeff_(move(coeff.scalarCoeff_)),
    tensorCoeff_(move(coeff.tensorCoeff_))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::LagrangianEqnBase& Foam::LagrangianSp<Foam::vector>::eqn() const
{
    return scalarCoeff_.eqn();
}


bool Foam::LagrangianSp<Foam::vector>::valid() const
{
    return
        scalarCoeff_.valid()
     || tensorCoeff_.valid();
}


void Foam::LagrangianSp<Foam::vector>::negate()
{
    scalarCoeff_.negate();
    tensorCoeff_.negate();
}


Foam::tmp<Foam::LagrangianCoeff<Foam::scalar, true>>
Foam::LagrangianSp<Foam::vector>::A() const
{
    return
        tensorCoeff_.valid()
      ? tmp<LagrangianCoeff<scalar, true>>
        (
            new LagrangianCoeff<scalar, true>
            (
                eqn(),
                (1.0/3.0)*tr(tensorCoeff_.S())
            )
        )
      : tmp<LagrangianCoeff<scalar, true>>(scalarCoeff_);
}


Foam::tmp<Foam::LagrangianCoeff<Foam::vector, false>>
Foam::LagrangianSp<Foam::vector>::H() const
{
    const LagrangianEqn<vector>& eqn =
        static_cast<const LagrangianEqn<vector>&>(this->eqn());

    return
        tmp<LagrangianCoeff<vector, false>>
        (
            tensorCoeff_.valid()
          ? new LagrangianCoeff<vector, false>
            (
                eqn,
                dev2(tensorCoeff_.S()) & eqn.psi()
            )
          : new LagrangianCoeff<vector, false>(eqn)
        );
}


Foam::tmp<Foam::LagrangianCoeff<Foam::vector, false>>
Foam::LagrangianSp<Foam::vector>::Su() const
{
    return Su(static_cast<const LagrangianEqn<vector>&>(this->eqn()));
}


Foam::tmp<Foam::LagrangianCoeff<Foam::vector, false>>
Foam::LagrangianSp<Foam::vector>::Su(const LagrangianEqn<vector>& eqn) const
{
    return
        tmp<LagrangianCoeff<vector, false>>
        (
            valid()
          ? new LagrangianCoeff<vector, false>
            (
                eqn,
                tensorCoeff_.valid()
              ? tensorCoeff_.S() & eqn.psi()
              : scalarCoeff_.S()*eqn.psi()
            )
          : new LagrangianCoeff<vector, false>(eqn)
        );
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const LagrangianCoeff<scalar, true>& coeff
)
{
    if (!coeff.valid()) return;
    operator+=(coeff.S());
}


void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const LagrangianCoeff<tensor, true>& coeff
)
{
    if (!coeff.valid()) return;
    operator+=(coeff.S());
}


void Foam::LagrangianSp<Foam::vector>::operator+=
(
    const LagrangianSp<vector>& Sp
)
{
    if (Sp.tensorCoeff_.valid())
    {
        operator+=(Sp.tensorCoeff_.S());
    }
    if (Sp.scalarCoeff_.valid())
    {
        operator+=(Sp.scalarCoeff_.S());
    }
}


void Foam::LagrangianSp<Foam::vector>::operator+=(const dimensioned<scalar>& dt)
{
    if (tensorCoeff_.valid())
    {
          tensorCoeff_ += tensor::I*dt;
    }
    else
    {
        scalarCoeff_ += dt;
    }
}


void Foam::LagrangianSp<Foam::vector>::operator+=(const dimensioned<tensor>& dt)
{
    makeTensor();
    tensorCoeff_ += dt;
}


void Foam::LagrangianSp<Foam::vector>::operator+=(const zero)
{}


void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const LagrangianCoeff<scalar, true>& coeff
)
{
    if (!coeff.valid()) return;
    operator-=(coeff.S());
}


void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const LagrangianCoeff<tensor, true>& coeff
)
{
    if (!coeff.valid()) return;
    operator-=(coeff.S());
}


void Foam::LagrangianSp<Foam::vector>::operator-=
(
    const LagrangianSp<vector>& Sp
)
{
    if (Sp.tensorCoeff_.valid())
    {
        operator-=(Sp.tensorCoeff_.S());
    }
    if (Sp.scalarCoeff_.valid())
    {
        operator-=(Sp.scalarCoeff_.S());
    }
}


void Foam::LagrangianSp<Foam::vector>::operator-=(const dimensioned<scalar>& dt)
{
    if (tensorCoeff_.valid())
    {
        tensorCoeff_ -= tensor::I*dt;
    }
    else
    {
        scalarCoeff_ -= dt;
    }
}


void Foam::LagrangianSp<Foam::vector>::operator-=(const dimensioned<tensor>& dt)
{
    makeTensor();
    tensorCoeff_ -= dt;
}


void Foam::LagrangianSp<Foam::vector>::operator-=(const zero)
{}


void Foam::LagrangianSp<Foam::vector>::operator*=(const dimensioned<scalar>& dt)
{
    scalarCoeff_ *= dt;
    tensorCoeff_ *= dt;
}


void Foam::LagrangianSp<Foam::vector>::operator/=(const dimensioned<scalar>& dt)
{
    scalarCoeff_ /= dt;
    tensorCoeff_ /= dt;
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubField<Foam::vector>> Foam::operator/
(
    const LagrangianCoeff<vector, false>& Su,
    const LagrangianSp<vector>& Sp
)
{
    if (Sp.tensorCoeff_.valid())
    {
        return inv(Sp.tensorCoeff_.S()) & Su.S();
    }
    else
    {
        return Su.S()/Sp.scalarCoeff_.S();
    }
}


// ************************************************************************* //
