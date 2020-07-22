/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "thermoParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::thermoParcelInjectionData::thermoParcelInjectionData(Istream& is)
:
    momentumParcelInjectionData(is)
{
    is.check("reading T");
    is >> T_;

    is.check("reading Cp");
    is >> Cp_;

    is.check("thermoParcelInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const thermoParcelInjectionData& data
)
{
    os << static_cast<const momentumParcelInjectionData&>(data);

    os << data.T_ << data.Cp_;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, thermoParcelInjectionData& data)
{
    is >> static_cast<momentumParcelInjectionData&>(data);

    is.check("reading T");
    is >> data.T_;

    is.check("reading Cp");
    is >> data.Cp_;

    is.check("operator(Istream&, thermoParcelInjectionData&)");

    return is;
}


// ************************************************************************* //
