/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

const Foam::NamedEnum<Foam::timeControl::timeControls, 9>
Foam::timeControl::timeControlNames_
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::timeControl::roundDown(const scalar t)
{
    return t < 0 ? label(t) - 1 : label(t);
}


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
    startTime_(time_.beginTime().value()),
    endTime_(vGreat),
    beginTime_(time_.beginTime().value()),
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
    dict.readIfPresent("startTime", time().userUnits(), startTime_);
    dict.readIfPresent("endTime", time().userUnits(), endTime_);
    dict.readIfPresent("beginTime", time().userUnits(), beginTime_);

    word controlName(prefix_ + "Control");
    word intervalName(prefix_ + "Interval");

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
        {
            interval_ = dict.lookup<scalar>(intervalName, time_.userUnits());
            break;
        }

        case timeControls::adjustableRunTime:
        {
            interval_ = dict.lookup<scalar>(intervalName, time_.userUnits());
            executionIndex_ =
                roundDown
                (
                    (
                        max(time_.value(), startTime_)
                      - beginTime_
                      - 0.5*time_.deltaTValue()
                    )
                   /interval_
                );
            break;
        }

        case timeControls::runTimes:
        {
            const word timesName(prefix_ + "Times");
            const word frequenciesName(prefix_ + "Frequencies");
            const bool repeat = dict.lookupOrDefault(prefix_ + "Repeat", false);

            timeDelta_ =
                dict.lookupOrDefault
                (
                    "timeDelta",
                    unitNone,
                    1e-3*time_.userDeltaTValue()
                );

            if (dict.found(timesName))
            {
                times_ = dict.lookup<scalarList>(timesName, unitNone);
            }
            else if (dict.found(frequenciesName))
            {
                List<Pair<scalar>> frequencies(dict.lookup(frequenciesName));

                const scalar userEndTime =
                    time_.timeToUserTime(time_.endTime().value());

                if (!repeat)
                {
                    frequencies.append
                    (
                        {userEndTime, frequencies.last().second()}
                    );
                }

                const scalar frequenciesDuration =
                    frequencies.last().first() - frequencies.first().first();

                DynamicList<scalar> times(1, frequencies[0].first());
                label i = 0;
                label repeati = 0;
                while (times[i] < userEndTime)
                {
                    for(label pi=0; pi<frequencies.size()-1; pi++)
                    {
                        while
                        (
                            times[i]
                          < frequencies[pi + 1].first()
                          + repeati*frequenciesDuration - timeDelta_
                        )
                        {
                            times(i + 1) = times[i] + frequencies[pi].second();
                            i++;
                        }
                    }
                    repeati++;
                }

                times_ = times;
            }
            else
            {
                FatalErrorInFunction
                    << "Undefined " << timesName
                    << " or " << frequenciesName << " for output control: "
                    << timeControlNames_[timeControl_] << nl
                    << exit(FatalError);
            }

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


bool Foam::timeControl::active() const
{
    return
        time_.value() >= startTime_ - 0.5*time_.deltaTValue()
     && time_.value() <= endTime_;
}


bool Foam::timeControl::execute()
{
    if (!active())
    {
        return false;
    }

    switch (timeControl_)
    {
        case timeControls::timeStep:
        {
            return
                (intervalSteps_ <= 1)
             || !(time_.timeIndex() % intervalSteps_);
            break;
        }

        case timeControls::writeTime:
        case timeControls::outputTime:
        {
            if
            (
                time_.writeTime()
             || time_.timeIndex() == time_.startTimeIndex()
            )
            {
                executionIndex_++;
                return !(executionIndex_ % intervalSteps_);
            }
            break;
        }

        case timeControls::runTime:
        case timeControls::adjustableRunTime:
        {
            const label executionIndex =
                roundDown
                (
                    (
                        time_.value()
                      - beginTime_
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
                << exit(FatalError);
            break;
        }
    }

    return false;
}


Foam::scalar Foam::timeControl::timeToNextAction()
{
    if (time_.value() > endTime_) return vGreat;

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
                (executionIndex_ + 1)*interval_
              - (time_.value() - beginTime_);
            break;
        }

        case timeControls::runTimes:
        {
            scalar realTimeToNextAction = vGreat;

            forAll(times_, i)
            {
                if
                (
                    time_.userTimeToTime(times_[i]) < startTime_
                 || time_.userTimeToTime(times_[i]) > endTime_
                ) continue;

                const scalar userTimeToThisAction =
                    times_[i] - time_.userTimeValue();

                if (userTimeToThisAction > timeDelta_)
                {
                    realTimeToNextAction =
                        min
                        (
                            realTimeToNextAction,
                            time_.userTimeToTime(userTimeToThisAction)
                        );
                }
            }

            return realTimeToNextAction;
            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Undefined output control: "
                << timeControlNames_[timeControl_] << nl
                << exit(FatalError);
            break;
        }
    }

    return vGreat;
}


// ************************************************************************* //
