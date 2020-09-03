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

#include "nonUniformTable1DThermophysicalFunction.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::thermophysicalFunctions::nonUniformTable1D::index
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
            << "    of nonUniformTable1D " << name_
            << exit(FatalError);
    }

    const scalar nd = (T - Tlow_)/deltaT_;
    const label j = nd;

    label i = jumpTable_[j];

    if (i < values_.size() - 1 && T > values_[i + 1].first())
    {
        i++;
    }

    return i;
}


// ************************************************************************* //
