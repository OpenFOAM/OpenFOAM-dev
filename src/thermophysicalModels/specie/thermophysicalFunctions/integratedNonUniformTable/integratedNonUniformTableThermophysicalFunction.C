/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "integratedNonUniformTableThermophysicalFunction.H"
#include "specie.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{
    defineTypeNameAndDebug(integratedNonUniformTable, 0);

    addToRunTimeSelectionTable
    (
        thermophysicalFunction,
        integratedNonUniformTable,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalFunctions::integratedNonUniformTable::
integratedNonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    nonUniformTable(name, dict),
    intf_(values().size()),
    intfByT_(values().size())
{
    intf_[0] = 0;
    intfByT_[0] = 0;

    for(label i = 0; i<intf_.size() - 1; i++)
    {
        intf_[i + 1] = intf_[i] + intfdT(0, values()[i + 1].first());
        intfByT_[i + 1] = intfByT_[i] + intfByTdT(0, values()[i + 1].first());
    }

    const scalar intfStd = intfdT(Pstd, Tstd);
    const scalar intfByTStd = intfByTdT(Pstd, Tstd);

    forAll(intf_, i)
    {
        intf_[i] -= intfStd;
        intfByT_[i] -= intfByTStd;
    }
}


Foam::thermophysicalFunctions::integratedNonUniformTable::
integratedNonUniformTable
(
    const dictionary& dict
)
:
    integratedNonUniformTable("values", dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::thermophysicalFunctions::integratedNonUniformTable::intfdT
(
    scalar p,
    scalar T
) const
{
    const label i = index(p, T);
    const scalar Ti = values()[i].first();
    const scalar fi = values()[i].second();
    const scalar dT = T - Ti;
    const scalar lambda = dT/(values()[i + 1].first() - Ti);

    return
        intf_[i]
      + (fi + 0.5*lambda*(values()[i + 1].second() - fi))*dT;
}


Foam::scalar Foam::thermophysicalFunctions::integratedNonUniformTable::intfByTdT
(
    scalar p,
    scalar T
) const
{
    const label i = index(p, T);
    const scalar Ti = values()[i].first();
    const scalar fi = values()[i].second();
    const scalar gradf =
        (values()[i + 1].second() - fi)/(values()[i + 1].first() - Ti);

    return
        intfByT_[i] + ((fi - gradf*Ti)*log(T/Ti) + gradf*(T - Ti));
}


void Foam::thermophysicalFunctions::integratedNonUniformTable::write
(
    Ostream& os
) const
{
    nonUniformTable::write(os);
}


// ************************************************************************* //
