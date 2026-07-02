/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "NSRDS.H"
#include "NSRDS5.H"
#include "dictionary.H"
#include "units.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::NSRDS<Specie>::NSRDS
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Specie(name, dict),
    a_(subDict.lookup<scalar>("a", units::none)),
    b_(subDict.lookup<scalar>("b", units::none)),
    c_(subDict.lookup<scalar>("c", dimTemperature)),
    d_(subDict.lookup<scalar>("d", dimless))
{
    if (subDict.found("type"))
    {
        const word functionType = subDict.lookup<word>("type");

        if (functionType != Function1s::NSRDS5::typeName)
        {
            FatalIOErrorInFunction(subDict)
                << "Density function for NSRDS equation of state "
                << "should be of type " << Function1s::NSRDS5::typeName
                << ". Function type given is " << functionType << "."
                << exit(FatalIOError);
        }
    }
}


template<class Specie>
Foam::NSRDS<Specie>::NSRDS(const word& name, const dictionary& dict)
:
    NSRDS(name, dict, dict.subDict("equationOfState"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::NSRDS<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    writeEntry
    (
        os,
        "equationOfState",
        dictionary::entries("a", a_, "b", b_, "c", c_, "d", d_)
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const NSRDS<Specie>& ns)
{
    ns.write(os);
    return os;
}


// ************************************************************************* //
