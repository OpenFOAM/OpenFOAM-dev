/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "OutputFilterFunctionObject.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::readControls()
{
    dict_.readIfPresent("timeStart", timeStart_);
    dict_.readIfPresent("timeEnd", timeEnd_);
    dict_.readIfPresent("nStepsToStartTimeChange", nStepsToStartTimeChange_);
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::active() const
{
    return
        time_.value() >= timeStart_
     && time_.value() <= timeEnd_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::OutputFilterFunctionObject<OutputFilter>::OutputFilterFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    dict_(dict),
    timeStart_(-VGREAT),
    timeEnd_(VGREAT),
    nStepsToStartTimeChange_
    (
        dict.lookupOrDefault("nStepsToStartTimeChange", 3)
    ),
    writeControl_(t, dict, "write"),
    evaluateControl_(t, dict, "evaluate"),
    filter_
    (
        name,
        time_.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        ),
        dict_
    )
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::execute
(
    const bool forceWrite
)
{
    if (active())
    {
        if (evaluateControl_.output())
        {
            filter_.execute();
        }

        if (forceWrite || writeControl_.output())
        {
            filter_.write();
        }
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::end()
{
    filter_.end();

    if (writeControl_.output())
    {
        filter_.write();
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::timeSet()
{
    if (active())
    {
        filter_.timeSet();
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::adjustTimeStep()
{
    if
    (
        active()
     && writeControl_.writeControl()
     == outputFilterControl::ocAdjustableRunTime
    )
    {
        const label  outputTimeIndex = writeControl_.outputTimeLastDump();
        const scalar writeInterval = writeControl_.writeInterval();

        scalar timeToNextWrite = max
        (
            0.0,
            (outputTimeIndex + 1)*writeInterval
          - (time_.value() - time_.startTime().value())
        );

        scalar deltaT = time_.deltaTValue();

        scalar nSteps = timeToNextWrite/deltaT - SMALL;

        // functionObjects modify deltaT within nStepsToStartTimeChange
        // NOTE: Potential problems arise if two function objects dump within
        // the same interval
        if (nSteps < nStepsToStartTimeChange_)
        {
            label nStepsToNextWrite = label(nSteps) + 1;

            scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

            // Adjust time step
            if (newDeltaT < deltaT)
            {
                deltaT = max(newDeltaT, 0.2*deltaT);
                const_cast<Time&>(time_).setDeltaT(deltaT, false);
            }
        }
    }

    return true;
}


template<class OutputFilter>
bool Foam::OutputFilterFunctionObject<OutputFilter>::read
(
    const dictionary& dict
)
{
    if (dict != dict_)
    {
        dict_ = dict;

        writeControl_.read(dict);
        evaluateControl_.read(dict);
        readControls();

        return true;
    }
    else
    {
        return false;
    }
}


template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::updateMesh
(
    const mapPolyMesh& mpm
)
{
    if (active())
    {
        filter_.updateMesh(mpm);
    }
}


template<class OutputFilter>
void Foam::OutputFilterFunctionObject<OutputFilter>::movePoints
(
    const polyMesh& mesh
)
{
    if (active())
    {
        filter_.movePoints(mesh);
    }
}


// ************************************************************************* //
