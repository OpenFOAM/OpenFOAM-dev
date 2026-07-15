/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "constTransport.H"
#include "dictionary.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::constTransport<Thermo>::constTransport
(
    const word& name,
    const dictionary& dict
)
:
    Thermo(name, dict)
{
    const dictionary& transportDict = dict.subDict("transport");

    mu_ = transportDict.lookup<scalar>("mu", dimensions::dynamicViscosity);

    const bool foundPr = transportDict.found("Pr");
    const bool foundKappa = transportDict.found("kappa");

    if (foundPr == foundKappa)
    {
        FatalIOErrorInFunction(dict)
            << "Either Pr or kappa must be specified, but not both."
            << exit(FatalIOError);
    }

    constPr_ = foundPr;
    rPr_ =
        constPr_
      ? 1/transportDict.lookup<scalar>("Pr", dimless)
      : NaN;
    kappa_ =
        constPr_
      ? NaN
      : transportDict.lookup<scalar>("kappa", dimensions::thermalConductivity);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::constTransport<Thermo>::constTransport::write(Ostream& os) const
{
    Thermo::write(os);

    const delimitDictionary delimit(os, "transport");
    writeEntry(os, "mu", mu_);
    if (constPr_)
    {
        writeEntry(os, "Pr", 1/rPr_);
    }
    else
    {
        writeEntry(os, "kappa", kappa_);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const constTransport<Thermo>& ct)
{
    ct.write(os);
    return os;
}


// ************************************************************************* //
