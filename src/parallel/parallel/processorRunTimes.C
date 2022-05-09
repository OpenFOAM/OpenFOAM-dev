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

#include "processorRunTimes.H"
#include "decompositionMethod.H"
#include "timeSelector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorRunTimes::processorRunTimes
(
    const word& name,
    const argList& args
)
:
    completeRunTime_(name, args),
    procRunTimes_
    (
        decompositionMethod::decomposeParDict(completeRunTime_)
       .lookup<int>("numberOfSubdomains")
    )
{
    forAll(procRunTimes_, proci)
    {
        procRunTimes_.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                completeRunTime_.rootPath(),
                completeRunTime_.caseName()
               /fileName(word("processor") + Foam::name(proci))
            )
        );

        procRunTimes_[proci].setTime(completeRunTime_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorRunTimes::~processorRunTimes()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::processorRunTimes::setTime
(
    const instant& inst,
    const label newIndex
)
{
    completeRunTime_.setTime(inst, newIndex);

    forAll(procRunTimes_, proci)
    {
        procRunTimes_[proci].setTime(inst, newIndex);
    }
}


Foam::instantList Foam::processorRunTimes::selectComplete(const argList& args)
{
    instantList timeDirs =
        timeSelector::selectIfPresent(completeRunTime_, args);

    forAll(procRunTimes_, proci)
    {
        procRunTimes_[proci].setTime(completeRunTime_);
    }

    return timeDirs;
}


Foam::instantList Foam::processorRunTimes::selectProc(const argList& args)
{
    instantList timeDirs =
        timeSelector::select0(procRunTimes_[0], args);

    completeRunTime_.setTime(procRunTimes_[0]);

    for (label proci = 1; proci < nProcs(); proci ++)
    {
        procRunTimes_[proci].setTime(procRunTimes_[0]);
    }

    return timeDirs;
}


// ************************************************************************* //
