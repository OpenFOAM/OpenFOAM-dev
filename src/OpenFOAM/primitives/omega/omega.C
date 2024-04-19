/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "omega.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omega::omega(const dictionary& dict)
:
    dimensionedScalar("omega", dimless/dimTime, NaN)
{
    const bool foundOmega = dict.found("omega");
    const bool foundRpm = dict.found("rpm");

    if (foundOmega && foundRpm)
    {
        FatalIOErrorInFunction(dict)
            << "Rotational speeds rpm and omega both defined in dictionary "
            << dict.name() << exit(FatalIOError);
    }

    if (!foundOmega && !foundRpm)
    {
        FatalIOErrorInFunction(dict)
            << "Neither rotational speed rpm or omega defined in dictionary "
            << dict.name() << exit(FatalIOError);
    }

    value() =
        foundOmega
      ? dict.lookup<scalar>("omega")
      : constant::mathematical::pi/30*dict.lookup<scalar>("rpm");
}


// ************************************************************************* //
