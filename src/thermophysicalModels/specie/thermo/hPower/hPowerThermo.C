/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "hPowerThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hPowerThermo<EquationOfState>::hPowerThermo(Istream& is)
:
    EquationOfState(is),
    n0_(readScalar(is)),
    Tref_(readScalar(is)),
    Hf_(readScalar(is))
{
    is.check("hPowerThermo::hPowerThermo(Istream& is)");
}


template<class EquationOfState>
Foam::hPowerThermo<EquationOfState>::hPowerThermo
(
    const dictionary& dict
)
:
    EquationOfState(dict),
    c0_(readScalar(dict.subDict("thermodynamics").lookup("C0"))),
    n0_(readScalar(dict.subDict("thermodynamics").lookup("n0"))),
    Tref_(readScalar(dict.subDict("thermodynamics").lookup("Tref"))),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf")))
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hPowerThermo<EquationOfState>& et
)
{
    os  << static_cast<const EquationOfState&>(et) << nl
        << "    " << et.c0_
        << tab << et.n0_
        << tab << et.Tref_
        << tab << et.Hf_;

    os << nl << "    ";

    os << endl;

    os.check
    (
        "operator<<(Ostream& os, const hPowerThermo<EquationOfState>& et)"
    );

    return os;
}


// ************************************************************************* //
