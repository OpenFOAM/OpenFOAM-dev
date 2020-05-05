/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "tableThermophysicalFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermophysicalFunctions
{
    defineTypeNameAndDebug(table, 0);
    addToRunTimeSelectionTable(thermophysicalFunction, table, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalFunctions::table::table(const dictionary& dict)
:
    dictName_(dict.name()),
    Tlow_(dict.lookup<scalar>("Tlow")),
    Thigh_(dict.lookup<scalar>("Thigh")),
    values_(dict.lookup("values"))
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
        deltaT_ = (Thigh_ - Tlow_)/(values_.size() - 1);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::thermophysicalFunctions::table::f(scalar p, scalar T) const
{
    const scalar nd = (T - Tlow_)/deltaT_;
    const label i = nd;

    if (nd < 0 || i > values_.size() - 2)
    {
        FatalErrorInFunction
            << "Temperature " << T << " out of range "
            << Tlow_ << " to " << Thigh_ << nl
            << "    of table " << dictName_
            << exit(FatalError);
    }

    const scalar Ti = Tlow_ + i*deltaT_;
    const scalar lambda = (T - Ti)/deltaT_;

    return values_[i] + lambda*(values_[i + 1] - values_[i]);
}


void Foam::thermophysicalFunctions::table::write(Ostream& os) const
{
    writeEntry(os, "Tlow", Tlow_);
    writeEntry(os, "Thigh", Thigh_);
    writeEntry(os, "values", values_);
}


// ************************************************************************* //
