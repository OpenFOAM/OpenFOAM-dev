/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "omega1.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::Function1<Foam::scalar>> Foam::Function1s::omega::init
(
    const Time& time,
    const dictionary& dict
)
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

    return
        Function1<scalar>::New
        (
            foundOmega ? "omega" : "rpm",
            tUnits_,
            foundOmega ? unitRadians/dimTime : units()["rpm"],
            dict
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::omega::omega(const Time& time, const dictionary& dict)
:
    tUnits_(time.userUnits()),
    omega_(init(time, dict))
{}


Foam::Function1s::omega::omega(const omega& o)
:
    tUnits_(o.tUnits_),
    omega_(o.omega_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function1s::omega::write(Ostream& os) const
{
    writeEntry
    (
        os,
        tUnits_,
        omega_->name() == "omega" ? unitRadians/dimTime : units()["rpm"],
        omega_()
    );
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

void Foam::Function1s::writeEntry(Ostream& os, const omega& o)
{
    o.write(os);
}


// ************************************************************************* //
