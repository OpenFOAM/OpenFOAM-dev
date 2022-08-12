/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

Foam::Function1s::omega::omega()
:
    rpm_(false),
    omegaFactor_(1)
{}


Foam::Function1s::omega::omega(const dictionary& dict)
:
    rpm_(dict.found("rpm")),
    omegaFactor_(rpm_ ? constant::mathematical::pi/30.0 : 1),
    omega_
    (
        rpm_
      ? Function1<scalar>::New("rpm", dict)
      : Function1<scalar>::New("omega", dict)
    )
{}


Foam::Function1s::omega::omega(const omega& o)
:
    rpm_(o.rpm_),
    omegaFactor_(o.omegaFactor_),
    omega_(o.omega_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::Function1s::omega::read(const dictionary& dict)
{
    rpm_ = dict.found("rpm");
    omegaFactor_ = rpm_ ? constant::mathematical::pi/30.0 : 1;

    omega_.reset
    (
        rpm_
      ? Function1<scalar>::New("rpm", dict).ptr()
      : Function1<scalar>::New("omega", dict).ptr()
    );

    return true;
}


void Foam::Function1s::omega::write(Ostream& os) const
{
    writeEntry(os, omega_());
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

void Foam::Function1s::writeEntry(Ostream& os, const omega& o)
{
    o.write(os);
}


// ************************************************************************* //
