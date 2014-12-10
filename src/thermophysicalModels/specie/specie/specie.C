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

#include "specie.H"
#include "IOstreams.H"
#include "constants.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

//- Universal gas constant (default in [J/(kmol K)])
const Foam::scalar Foam::specie::RR = constant::physicoChemical::R.value()*1000;

//- Standard pressure (default in [Pa])
const Foam::scalar Foam::specie::Pstd = constant::standard::Pstd.value();

//- Standard temperature (default in [K])
const Foam::scalar Foam::specie::Tstd = constant::standard::Tstd.value();

namespace Foam
{
defineTypeNameAndDebug(specie, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(Istream& is)
:
    name_(is),
    nMoles_(readScalar(is)),
    molWeight_(readScalar(is))
{
    is.check("specie::specie(Istream& is)");
}


Foam::specie::specie(const dictionary& dict)
:
    name_(dict.dictName()),
    nMoles_(readScalar(dict.subDict("specie").lookup("nMoles"))),
    molWeight_(readScalar(dict.subDict("specie").lookup("molWeight")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specie::write(Ostream& os) const
{
    dictionary dict("specie");
    dict.add("nMoles", nMoles_);
    dict.add("molWeight", molWeight_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    os  << st.name_ << tab
        << st.nMoles_ << tab
        << st.molWeight_;

    os.check("Ostream& operator<<(Ostream& os, const specie& st)");
    return os;
}


// ************************************************************************* //
