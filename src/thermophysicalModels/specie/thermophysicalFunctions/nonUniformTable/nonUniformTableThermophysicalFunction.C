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

#include "nonUniformTableThermophysicalFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{
    defineTypeNameAndDebug(nonUniformTable, 0);

    addToRunTimeSelectionTable
    (
        thermophysicalFunction,
        nonUniformTable,
        dictionary
    );
}
}

template<>
const char* const Foam::Tuple2<Foam::scalar, Foam::scalar>::typeName
(
    "Tuple2<scalar,scalar>"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalFunctions::nonUniformTable::nonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    Tlow_(great),
    Thigh_(-great),
    values_(dict.lookup(name)),
    deltaT_(great)
{
    if (values_.size() < 2)
    {
        FatalIOErrorInFunction(dict)
            << "Table " << nl
            << "    " << name_ << nl
            << "    has less than 2 entries."
            << exit(FatalIOError);
    }
    else
    {
        Tlow_ = values_.first().first();
        Thigh_ = values_.last().first();

        for(label i = 1; i<values_.size(); i++)
        {
            deltaT_ = min(deltaT_, values_[i].first() - values_[i - 1].first());
        }

        deltaT_ *= 0.9;

        jumpTable_.setSize((Thigh_ - Tlow_)/deltaT_ + 1);

        label i = 0;
        forAll(jumpTable_, j)
        {
            const scalar T = Tlow_ + j*deltaT_;

            if (T > values_[i + 1].first())
            {
                i++;
            }

            jumpTable_[j] = i;
        }
    }
}


Foam::thermophysicalFunctions::nonUniformTable::nonUniformTable
(
    const dictionary& dict
)
:
    nonUniformTable("values", dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::thermophysicalFunctions::nonUniformTable::f
(
    scalar p,
    scalar T
) const
{
    const label i = index(p, T);
    const scalar Ti = values_[i].first();
    const scalar lambda = (T - Ti)/(values_[i + 1].first() - Ti);

    return
        values_[i].second()
      + lambda*(values_[i + 1].second() - values_[i].second());
}


Foam::scalar Foam::thermophysicalFunctions::nonUniformTable::dfdT
(
    scalar p,
    scalar T
) const
{
    const label i = index(p, T);

    return
        (values_[i + 1].second() - values_[i].second())
       /(values_[i + 1].first() - values_[i].first());
}


void Foam::thermophysicalFunctions::nonUniformTable::write(Ostream& os) const
{
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
