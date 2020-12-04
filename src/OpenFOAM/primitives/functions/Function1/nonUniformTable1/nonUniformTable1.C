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

#include "nonUniformTable1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    makeScalarFunction1(nonUniformTable)
}
}

template<>
const char* const Foam::Tuple2<Foam::scalar, Foam::scalar>::typeName
(
    "Tuple2<scalar,scalar>"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::nonUniformTable::nonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction1<scalar, nonUniformTable>(name),
    low_(great),
    high_(-great),
    values_(dict.lookup("values")),
    delta_(great)
{
    if (values_.size() < 2)
    {
        FatalIOErrorInFunction(dict)
            << "Table " << nl
            << "    " << name << nl
            << "    has less than 2 entries."
            << exit(FatalIOError);
    }
    else
    {
        low_ = values_.first().first();
        high_ = values_.last().first();

        for(label i = 1; i<values_.size(); i++)
        {
            delta_ = min(delta_, values_[i].first() - values_[i - 1].first());
        }

        delta_ *= 0.9;

        jumpTable_.setSize((high_ - low_)/delta_ + 1);

        label i = 0;
        forAll(jumpTable_, j)
        {
            const scalar x = low_ + j*delta_;

            if (x > values_[i + 1].first())
            {
                i++;
            }

            jumpTable_[j] = i;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1s::nonUniformTable::value
(
    scalar x
) const
{
    const label i = index(x);
    const scalar xi = values_[i].first();
    const scalar lambda = (x - xi)/(values_[i + 1].first() - xi);

    return
        values_[i].second()
      + lambda*(values_[i + 1].second() - values_[i].second());
}


Foam::scalar Foam::Function1s::nonUniformTable::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return 0;
}


Foam::scalar Foam::Function1s::nonUniformTable::dfdT
(
    scalar T
) const
{
    const label i = index(T);

    return
        (values_[i + 1].second() - values_[i].second())
       /(values_[i + 1].first() - values_[i].first());
}


void Foam::Function1s::nonUniformTable::write(Ostream& os) const
{
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
