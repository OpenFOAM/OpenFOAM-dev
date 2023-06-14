/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "timeControl.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<timeControl::timeControls, 9>::
    names[] =
    {
        "timeStep",
        "writeTime",
        "outputTime",
        "adjustableRunTime",
        "runTime",
        "runTimes",
        "clockTime",
        "cpuTime",
        "none"
    };
}

const Foam::NamedEnum<Foam::timeControl::timeControls, 9>
    Foam::timeControl::timeControlNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeControl::timeControl
(
    const Time& t,
    const dictionary& dict,
    const word& prefix
)
:
    time_(t),
    prefix_(prefix),
    timeControl_(timeControls::timeStep),
    intervalSteps_(0),
    interval_(-1),
    timeDelta_(0),
    executionIndex_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeControl::~timeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeControl::read(const dictionary& dict)
{
    word controlName(prefix_ + "Control");
    word intervalName(prefix_ + "Interval");
    const word timesName(prefix_ + "Times");

    // For backward compatibility support the deprecated 'outputControl' option
    // now superseded by 'writeControl' for compatibility with Time
    if (prefix_ == "write" && dict.found("outputControl"))
    {
        IOWarningInFunction(dict)
            << "Using deprecated 'outputControl'" << nl
            << "    Please use 'writeControl' with 'writeInterval'"
            << endl;

        // Change to the old names for this option
        controlName = "outputControl";
        intervalName = "outputInterval";
    }

    if (dict.found(controlName))
    {
        timeControl_ = timeControlNames_.read(dict.lookup(controlName));
    }
    else
    {
        timeControl_ = timeControls::timeStep;
    }

    switch (timeControl_)
    {
        case timeControls::timeStep:
        {
            intervalSteps_ = dict.lookupOrDefault<label>(intervalName, 0);
            break;
        }

        case timeControls::writeTime:
        case timeControls::outputTime:
        {
            intervalSteps_ = dict.lookupOrDefault<label>(intervalName, 1);
            break;
        }

        case timeControls::clockTime:
        case timeControls::runTime:
        case timeControls::cpuTime:
        case timeControls::adjustableRunTime:
        {
            interval_ = time_.userTimeToTime(dict.lookup<scalar>(intervalName));

            if (timeControl_ == timeControls::adjustableRunTime)
            {
                executionIndex_ = label
                (
                    (
                        (time_.value() - time_.beginTime().value())
                      + 0.5*time_.deltaTValue()
                    )
                    /interval_
                );
            }

            break;
        }

        case timeControls::runTimes:
        {
            times_ = dict.lookup<scalarList>(timesName);
            timeDelta_ = dict.lookupOrDefault("timeDelta", 1e-6);

            forAll(times_, i)
            {
                timeIndices_.insert
                (
                    int64_t((times_[i] + timeDelta_/2.0)/timeDelta_)
                );
            }

            intervalSteps_ = dict.lookupOrDefault<label>(intervalName, 1);
            break;
        }

        default:
        {
            break;
        }
    }
}


bool Foam::timeControl::execute()
{
    switch (timeControl_)
    {
        case timeControls::timeStep:
        {
            return
            (
                (intervalSteps_ <= 1)
             || !(time_.timeIndex() % intervalSteps_)
            );
            break;
        }

        case timeControls::writeTime:
        case timeControls::outputTime:
        {
            if (time_.writeTime())
            {
                executionIndex_++;
                return !(executionIndex_ % intervalSteps_);
            }
            break;
        }

        case timeControls::runTime:
        case timeControls::adjustableRunTime:
        {
            const label executionIndex = label
            (
                (
                    (time_.value() - time_.beginTime().value())
                  + 0.5*time_.deltaTValue()
                )
               /interval_
            );

            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case timeControls::runTimes:
        {
            return timeIndices_.found
            (
                (time_.userTimeValue() + timeDelta_/2)/timeDelta_
            );

            break;
        }

        case timeControls::cpuTime:
        {
            const label executionIndex = label
            (
                returnReduce(time_.elapsedCpuTime(), maxOp<double>())
               /interval_
            );
            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case timeControls::clockTime:
        {
            const label executionIndex = label
            (
                returnReduce(label(time_.elapsedClockTime()), maxOp<label>())
               /interval_
            );
            if (executionIndex > executionIndex_)
            {
                executionIndex_ = executionIndex;
                return true;
            }
            break;
        }

        case timeControls::none:
        {
            return false;
        }

        default:
        {
            FatalErrorInFunction
                << "Undefined output control: "
                << timeControlNames_[timeControl_] << nl
                << abort(FatalError);
            break;
        }
    }

    return false;
}


Foam::scalar Foam::timeControl::timeToNextAction()
{
    switch (timeControl_)
    {
        case timeControls::timeStep:
        case timeControls::writeTime:
        case timeControls::outputTime:
        case timeControls::runTime:
        case timeControls::cpuTime:
        case timeControls::clockTime:
        case timeControls::none:
        {
            return vGreat;
            break;
        }

        case timeControls::adjustableRunTime:
        {
            return
                max
                (
                    0.0,
                    (executionIndex_ + 1)*interval_
                  - (time_.value() - time_.beginTime().value())
                );
            break;
        }

        case timeControls::runTimes:
        {
            if (time_.userTimeValue() + timeDelta_ < times_.last())
            {
                forAll(times_, i)
                {
                    if (times_[i] > time_.userTimeValue() + timeDelta_)
                    {
                        return time_.userTimeToTime
                        (
                            times_[i] - time_.userTimeValue()
                        );
                    }
                }
            }

            return vGreat;

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Undefined output control: "
                << timeControlNames_[timeControl_] << nl
                << abort(FatalError);
            break;
        }
    }

    return vGreat;
}


// ************************************************************************* //
