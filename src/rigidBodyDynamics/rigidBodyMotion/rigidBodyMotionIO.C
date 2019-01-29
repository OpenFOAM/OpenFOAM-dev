/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "rigidBodyMotion.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::RBD::rigidBodyMotion::read(const dictionary& dict)
{
    rigidBodyModel::read(dict);

    aRelax_ = dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0);
    aDamp_ = dict.lookupOrDefault<scalar>("accelerationDamping", 1.0);
    report_ = dict.lookupOrDefault<Switch>("report", false);

    return true;
}


void Foam::RBD::rigidBodyMotion::write(Ostream& os) const
{
    rigidBodyModel::write(os);

    writeEntry(os, "accelerationRelaxation", aRelax_);
    writeEntry(os, "accelerationDamping", aDamp_);
    writeEntry(os, "report", report_);
}


// ************************************************************************* //
