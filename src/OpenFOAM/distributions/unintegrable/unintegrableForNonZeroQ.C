/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "unintegrableForNonZeroQ.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::unintegrableForNonZeroQ::Phi
(
    const label q,
    const scalarField& x
) const
{
    if (q == 0)
    {
        return PhiForZeroQ(x);
    }
    else
    {
        return unintegrable::Phi(q, x);
    }
}


Foam::Pair<Foam::scalar> Foam::distributions::unintegrableForNonZeroQ::Phi01
(
    const label q
) const
{
    if (q == 0)
    {
        const scalarField x(scalarList({min(), max()}));
        const scalarField Phi(this->PhiForZeroQ(x));
        return Pair<scalar>(Phi.first(), Phi.last());
    }
    else
    {
        return unintegrable::Phi01(q);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::unintegrableForNonZeroQ::sample() const
{
    if (q() == 0)
    {
        return sampleForZeroQ();
    }
    else
    {
        return unintegrable::sample();
    }
}


Foam::tmp<Foam::scalarField>
Foam::distributions::unintegrableForNonZeroQ::CDF(const scalarField& x) const
{
    if (q() == 0)
    {
        const Pair<scalar>& Phi01 = this->Phi01();
        const scalarField xClip(Foam::min(Foam::max(x, min()), max()));
        return (PhiForZeroQ(xClip) - Phi01[0])/(Phi01[1] - Phi01[0]);
    }
    else
    {
        return unintegrable::CDF(x);
    }
}


// ************************************************************************* //
