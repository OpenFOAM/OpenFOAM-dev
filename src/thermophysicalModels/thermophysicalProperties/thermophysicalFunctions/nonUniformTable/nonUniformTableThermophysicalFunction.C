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
    const dictionary& dict
)
:
    dictName_(dict.name()),
    Tlow_(great),
    Thigh_(-great),
    values_(dict.lookup("values")),
    deltaT_(great)
{
    if (values_.size() < 2)
    {
        FatalErrorInFunction
            << "Table " << nl
            << "    " << dictName_ << nl
            << "    has less than 2 entries."
            << exit(FatalError);
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::thermophysicalFunctions::nonUniformTable::f
(
    scalar p,
    scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        FatalErrorInFunction
            << "Temperature " << T << " out of range "
            << Tlow_ << " to " << Thigh_ << nl
            << "    of nonUniformTable " << dictName_
            << exit(FatalError);
    }

    const scalar nd = (T - Tlow_)/deltaT_;
    const label j = nd;
    const label i = jumpTable_[j];

    const scalar Ti = values_[i].first();
    const scalar lambda = (T - Ti)/(values_[i + 1].first() - Ti);

    return
        values_[i].second()
      + lambda*(values_[i + 1].second() - values_[i].second());
}


void Foam::thermophysicalFunctions::nonUniformTable::write(Ostream& os) const
{
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
