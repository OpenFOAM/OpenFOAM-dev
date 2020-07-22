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

#include "momentumParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::momentumParcelInjectionData::momentumParcelInjectionData(Istream& is)
{
    is.check("reading (Px Py Pz)");
    is >> x_;

    is.check("reading (Ux Uy Uz)");
    is >> U_;

    is.check("reading d");
    is >> d_;

    is.check("reading rho");
    is >> rho_;

    is.check("reading mDot");
    is >> mDot_;

    is.check("momentumParcelInjectionData(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const momentumParcelInjectionData& data
)
{
    os << data.x_ << data.U_ << data.d_ << data.rho_ << data.mDot_;

    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, momentumParcelInjectionData& data)
{
    is.check("reading (Px Py Pz)");
    is >> data.x_;

    is.check("reading (Ux Uy Uz)");
    is >> data.U_;

    is.check("reading d");
    is >> data.d_;

    is.check("reading rho");
    is >> data.rho_;

    is.check("reading mDot");
    is >> data.mDot_;

    is.check("operator(Istream&, momentumParcelInjectionData&)");

    return is;
}


// ************************************************************************* //
