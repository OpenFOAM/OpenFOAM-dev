/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "eConstThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::eConstThermo<EquationOfState>::eConstThermo(Istream& is)
:
    EquationOfState(is),
    Cv_(readScalar(is)),
    Hf_(readScalar(is))
{
    is.check("eConstThermo<EquationOfState>::eConstThermo(Istream&)");

    Cv_ *= this->W();
    Hf_ *= this->W();
}


template<class EquationOfState>
Foam::eConstThermo<EquationOfState>::eConstThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Cv_(readScalar(dict.subDict("thermodynamics").lookup("Cv"))),
    Hf_(readScalar(dict.subDict("thermodynamics").lookup("Hf")))
{
    Cv_ *= this->W();
    Hf_ *= this->W();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::eConstThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Cv", Cv_/this->W());
    dict.add("Hf", Hf_/this->W());
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const eConstThermo<EquationOfState>& ct
)
{
    os  << static_cast<const EquationOfState&>(ct) << tab
        << ct.Cv_/ct.W() << tab << ct.Hf_/ct.W();

    os.check("Ostream& operator<<(Ostream&, const eConstThermo&)");
    return os;
}


// ************************************************************************* //
