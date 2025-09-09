/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::PtrList<Foam::Time>
Foam::processorRunTimes::processorRunTimes::procRunTimes
(
    const Time& completeRunTime,
    const nProcsFrom npf
)
{
    label nProcs;
    switch (npf)
    {
        case nProcsFrom::decomposeParDict:
            nProcs =
                decompositionMethod::decomposeParDict(completeRunTime)
               .lookup<int>("numberOfSubdomains");
            break;

        case nProcsFrom::fileHandler:
            nProcs =
                fileHandler().nProcs(completeRunTime.path());
            break;
    }

    PtrList<Time> result(nProcs);

    forAll(result, proci)
    {
        result.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                completeRunTime.rootPath(),
                completeRunTime.caseName()
               /fileName(word("processor") + Foam::name(proci))
            )
        );

        result[proci].setTime(completeRunTime);
    }

    return result;
}


Foam::autoPtr<Foam::Time> Foam::processorRunTimes::proc0RunTime
(
    const Time& completeRunTime,
    const PtrList<Time>& processorRunTimes
)
{
    return
        autoPtr<Time>
        (
            processorRunTimes.empty()
          ? new Time
            (
                Time::controlDictName,
                completeRunTime.rootPath(),
                completeRunTime.caseName()
               /fileName(word("processor") + Foam::name(label(0)))
            )
          : nullptr
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorRunTimes::processorRunTimes
(
    const word& name,
    const argList& args,
    const bool enableFunctionObjects,
    const nProcsFrom npf
)
:
    completeRunTime_(name, args, enableFunctionObjects),
    procRunTimes_(procRunTimes(completeRunTime_, npf)),
    proc0RunTime_(proc0RunTime(completeRunTime_, procRunTimes_))
{}


Foam::processorRunTimes::processorRunTimes
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const bool enableFunctionObjects,
    const nProcsFrom npf
)
:
    completeRunTime_(name, rootPath, caseName, enableFunctionObjects),
    procRunTimes_(procRunTimes(completeRunTime_, npf)),
    proc0RunTime_(proc0RunTime(completeRunTime_, procRunTimes_))
{}


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

    proc0TimeRef().setTime(inst, newIndex);

    for (label proci = 1; proci < nProcs(); proci ++)
    {
        procRunTimes_[proci].setTime(inst, newIndex);
    }
}


Foam::instantList Foam::processorRunTimes::selectComplete(const argList& args)
{
    const instantList timeDirs =
        timeSelector::selectIfPresent(completeRunTime_, args);

    proc0TimeRef().setTime(completeRunTime_);

    for (label proci = 1; proci < nProcs(); proci ++)
    {
        procRunTimes_[proci].setTime(completeRunTime_);
    }

    return timeDirs;
}


Foam::instantList Foam::processorRunTimes::selectProc(const argList& args)
{
    const instantList timeDirs =
        timeSelector::select0(proc0TimeRef(), args);

    completeRunTime_.setTime(proc0TimeRef());

    for (label proci = 1; proci < nProcs(); proci ++)
    {
        procRunTimes_[proci].setTime(proc0TimeRef());
    }

    return timeDirs;
}


// ************************************************************************* //
