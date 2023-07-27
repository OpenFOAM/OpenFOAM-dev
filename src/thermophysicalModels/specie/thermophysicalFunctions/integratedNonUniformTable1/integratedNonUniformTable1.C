/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "integratedNonUniformTable1.H"
#include "specie.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addScalarFunction1(integratedNonUniformTable);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::integratedNonUniformTable::integratedNonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    NonUniformTable<scalar>(name, dict),
    intf_(values().size()),
    intfByT_(values().size())
{
    intf_[0] = 0;
    intfByT_[0] = 0;

    for(label i = 0; i<intf_.size() - 1; i++)
    {
        intf_[i + 1] = intfdT(values()[i + 1].first());
        intfByT_[i + 1] = intfByTdT(values()[i + 1].first());
    }

    const scalar intfStd = intfdT(Tstd);
    const scalar intfByTStd = intfByTdT(Tstd);

    forAll(intf_, i)
    {
        intf_[i] -= intfStd;
        intfByT_[i] -= intfByTStd;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1s::integratedNonUniformTable::intfdT
(
    scalar T
) const
{
    const label i = index(T);
    const scalar Ti = values()[i].first();
    const scalar fi = values()[i].second();
    const scalar dT = T - Ti;
    const scalar lambda = dT/(values()[i + 1].first() - Ti);

    return
        intf_[i]
      + (fi + 0.5*lambda*(values()[i + 1].second() - fi))*dT;
}


Foam::scalar Foam::Function1s::integratedNonUniformTable::intfByTdT
(
    scalar T
) const
{
    const label i = index(T);
    const scalar Ti = values()[i].first();
    const scalar fi = values()[i].second();
    const scalar gradf =
        (values()[i + 1].second() - fi)/(values()[i + 1].first() - Ti);

    return
        intfByT_[i] + ((fi - gradf*Ti)*log(T/Ti) + gradf*(T - Ti));
}


void Foam::Function1s::integratedNonUniformTable::write
(
    Ostream& os
) const
{
    NonUniformTable<scalar>::write(os);
}


// ************************************************************************* //
