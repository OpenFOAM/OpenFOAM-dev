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

#include "data.H"
#include "Time.H"
#include "solverPerformance.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::data::debug(Foam::debug::debugSwitch("data", false));

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::data::data(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "data",
            obr.time().system(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    prevTimeIndex_(0)
{
    set("solverPerformance", dictionary());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::data::solverPerformanceDict() const
{
    return subDict("solverPerformance");
}


void Foam::data::setSolverPerformance
(
    const word& name,
    const solverPerformance& sp
) const
{
    dictionary& dict = const_cast<dictionary&>(solverPerformanceDict());

    List<solverPerformance> perfs;

    if (prevTimeIndex_ != this->time().timeIndex())
    {
        // reset solver performance between iterations
        prevTimeIndex_ = this->time().timeIndex();
        dict.clear();
    }
    else
    {
        dict.readIfPresent(name, perfs);
    }

    // append to list
    perfs.setSize(perfs.size()+1, sp);

    dict.set(name, perfs);
}


void Foam::data::setSolverPerformance
(
    const solverPerformance& sp
) const
{
    setSolverPerformance(sp.fieldName(), sp);
}


// ************************************************************************* //
