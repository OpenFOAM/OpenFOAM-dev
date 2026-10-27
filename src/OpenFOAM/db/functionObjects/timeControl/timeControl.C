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

#include "timeControl.H"
#include "Function1UnitSets.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Create a unit conversion object for tuples of scalar intervals
Function1s::unitSets tupleUnits(const unitSet& x, const unitSet& value)
{
    return Function1s::unitSets({x, value});
}


//- Create a unit conversion object for tuples of integer intervals
const unitSet& tupleUnits(const unitSet& x, const nil&)
{
    return x;
}


//- Read a tuple of scalar intervals
template<class Convert>
void convert
(
    Tuple2<scalar, token>& t,
    const Function1s::unitSets& units,
    const Convert& convert
)
{
    convert(t.first(), units.x);

    if (t.second().isNumber())
    {
        scalar s = t.second().number();
        convert(s, units.value);
        t.second() = token(s, t.second().lineNumber());
    }
}


//- Read a tuple of integer intervals
template<class Convert>
void convert
(
    Tuple2<scalar, token>& t,
    const unitSet& units,
    const Convert& convert
)
{
    convert(t.first(), units);
}


//- Get the value from a token
template<class Type>
Type tokenValue(const token& t);


//- Get the scalar value from a token
template<>
scalar tokenValue<scalar>(const token& t)
{
    return t.number();
}


//- Get the integer value from a token
template<>
label tokenValue<label>(const token& t)
{
    return t.labelToken();
}

}


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::NamedEnum<Foam::timeControl::timeControls, 8>
Foam::timeControl::timeControlNames_
{
    "timeStep",
    "writeTime",
    "outputTime",
    "adjustableRunTime",
    "runTime",
    "clockTime",
    "cpuTime",
    "none"
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class Units>
void Foam::timeControl::readIntervals
(
    const word& intervalName,
    const dictionary& dict,
    const Units& units,
    List<Type>& intervalSizes
)
{
    const bool haveInterval = dict.found(intervalName);
    const bool haveIntervals = dict.found(intervalName + "s");

    if (haveInterval == haveIntervals)
    {
        FatalIOErrorInFunction(dict)
            << (haveInterval ? "both" : "neither") << " of keywords "
            << intervalName << " " << (haveInterval ? "and" : "or") << " "
            << intervalName << "s defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    if (haveInterval)
    {
        intervalTimes_.resize(1);
        intervalTimes_[0] = startTime_;

        intervalSizes.resize(1);
        intervalSizes[0] = dict.lookup<Type>(intervalName, units);

        repeat_ = false;
    }
    else
    {
        List<Tuple2<scalar, token>> intervalTimeAndTokens =
            dict.lookup<List<Tuple2<scalar, token>>>
            (
                intervalName + "s",
                tupleUnits(time_.userUnits(), units)
            );

        intervalTimes_.resize(intervalTimeAndTokens.size());
        intervalSizes.resize(intervalTimeAndTokens.size());

        forAll(intervalTimeAndTokens, intervali)
        {
            const Tuple2<scalar, token>& timeAndToken =
                intervalTimeAndTokens[intervali];

            intervalTimes_[intervali] = timeAndToken.first();

            // Check monotonicity
            if
            (
                intervali
             && timeAndToken.first()
             <= intervalTimeAndTokens[intervali - 1].first()
            )
            {
                FatalIOErrorInFunction(dict)
                    << intervalName + "s do not increase monotonically"
                    << exit(FatalIOError);
            }

            auto error = [&]()
            {
                FatalIOErrorInFunction(dict)
                    << intervalName << " '" << timeAndToken.second()
                    << "' not recognised. Valid non-numeric values are "
                    << "'none', 'once' and 'repeat'." << exit(FatalIOError);
            };

            // Numeric interval
            if (timeAndToken.second().isNumber())
            {
                intervalSizes[intervali] =
                    tokenValue<Type>(timeAndToken.second());

                if (intervalSizes[intervali] <= 0)
                {
                    FatalIOErrorInFunction(dict)
                        << intervalName << " '" << intervalSizes[intervali]
                        << "' is not positive" << exit(FatalIOError);
                }
            }

            // Keyword interval
            else if (timeAndToken.second().isWord())
            {
                const word& key = timeAndToken.second().wordToken();

                if (key == "none")
                {
                    intervalSizes[intervali] = pTraits<Type>::max;
                }
                else if (key == "once")
                {
                    intervalSizes[intervali] = pTraits<Type>::min;
                }
                else if (key == "repeat")
                {
                    if (intervali != intervalTimeAndTokens.size() - 1)
                    {
                        FatalIOErrorInFunction(dict)
                            << intervalName << " '" << key
                            << "' must come last" << exit(FatalIOError);
                    }

                    intervalSizes[intervali] = 0;
                    repeat_ = true;
                }
                else
                {
                    error();
                }
            }

            // Unrecognised interval
            else
            {
                error();
            }
        }
    }
}


template<class Type, class Units>
void Foam::timeControl::readIntervals
(
    const word& intervalName,
    const dictionary& dict,
    const Units& units,
    const bool allowDefault,
    const Type& defaultInterval,
    List<Type>& intervalSizes
)
{
    const bool haveInterval = dict.found(intervalName);
    const bool haveIntervals = dict.found(intervalName + "s");

    if (allowDefault && !haveInterval && !haveIntervals)
    {
        intervalTimes_.resize(1);
        intervalTimes_[0] = startTime_;

        intervalSizes.resize(1);
        intervalSizes[0] = defaultInterval;

        repeat_ = false;
    }
    else
    {
        readIntervals<Type>
        (
            intervalName,
            dict,
            units,
            intervalSizes
        );
    }
}


Foam::scalar Foam::timeControl::repeating(const scalar t) const
{
    if (repeat_)
    {
        const scalar period =
            intervalTimes_.last() - intervalTimes_.first();

        return label(max(t - intervalTimes_.first(), 0)/period)*period;
    }
    else
    {
        return 0;
    }
}


Foam::scalar Foam::timeControl::after(const scalar t) const
{
    return t == 0 ? rootVSmall : t + rootSmall*mag(t);
}


template<class Type>
void Foam::timeControl::update
(
    const List<Type>& intervalSizes,
    const scalar t,
    const Type& x0,
    const Type& x,
    label& intervalIndex,
    label& executionIndex,
    Type& offset
) const
{
    // Find the current interval
    const label intervalIndex0 = intervalIndex;
    intervalIndex = max(intervalIndex, 0);
    const scalar t0 = repeating(t);
    if (repeat_)
    {
        for
        (
            ;
            max(t - t0, 0) < intervalTimes_[intervalIndex]
         || max(t - t0, 0) >= intervalTimes_[intervalIndex + 1];
            intervalIndex = (intervalIndex + 1) % (intervalTimes_.size() - 1)
        );
    }
    else
    {
        for
        (
            ;
            intervalIndex < intervalTimes_.size() - 1
         && t >= intervalTimes_[intervalIndex + 1];
         ++ intervalIndex
        );
    }

    // If in a new interval then compute the new execution index and offset
    if (intervalIndex != intervalIndex0)
    {
        const Type size = intervalSizes[intervalIndex];

        if (size == pTraits<Type>::max)
        {
            // Don't write
            executionIndex = 0;
            offset = 0;
        }
        else if (size == pTraits<Type>::min)
        {
            // Write once
            executionIndex = -1;
            offset = 0;
        }
        else
        {
            // Set the execution index and offset. Note: Round down, rather
            // than rounding towards zero, as is standard for integer division.
            const label executionIndex0 = executionIndex;
            const Type offset0 = offset;
            executionIndex = label((x - x0)/size);
            if (x - x0 < 0) executionIndex --;
            offset = 0;

            // If we have two valid/positive interval sizes then set the offset
            // to ensure a gap bounded by the sizes
            if (intervalIndex0 >= 0 && intervalSizes[intervalIndex0] > 0)
            {
                const Type size0 = intervalSizes[intervalIndex0];

                const Type xDelta =
                    ((executionIndex + 1)*size + x0 + offset)
                  - (executionIndex0*size0 + x0 + offset0);

                const Type minSize = min(size0, size);
                const Type maxSize = max(size0, size);

                if (xDelta < minSize) offset = minSize - xDelta;
                if (xDelta > maxSize) offset = maxSize - xDelta;
            }
        }
    }
}


template<class Type>
bool Foam::timeControl::execute
(
    const List<Type>& intervalSizes,
    const scalar t,
    const Type& x0,
    const Type& x,
    const Type& dx,
    Type& offset
)
{
    // Not yet in range
    if (t < intervalTimes_.first()) return false;

    // Update the interval
    update
    (
        intervalSizes,
        t,
        x0,
        x,
        intervalIndex_,
        executionIndex_,
        offset
    );

    // Determine the current execution index. Note: Round down, rather
    // than rounding towards zero, as is standard for integer division.
    const label executionIndex0 = executionIndex_;
    const Type size = mag(intervalSizes[intervalIndex_]);
    executionIndex_ = label((x + dx - x0 - offset)/size);
    if (x + dx - x0 - offset < 0) executionIndex_ --;

    // Return whether or not the index changed
    return executionIndex_ > executionIndex0;
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
    writeTimeIndex_(0),
    intervalTimes_(),
    intervalCounts_(),
    intervalPeriods_(),
    repeat_(false),
    intervalIndex_(labelMin),
    executionIndex_(labelMin),
    offsetCount_(0),
    offsetPeriod_(0)
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

    timeControl_ =
        timeControlNames_.lookupOrDefault
        (
            controlName,
            dict,
            timeControls::timeStep
        );

    dict.readIfPresent("startTime", time().userUnits(), startTime_);
    dict.readIfPresent("endTime", time().userUnits(), endTime_);
    dict.readIfPresent("beginTime", time().userUnits(), beginTime_);

    writeTimeIndex_ = 0;

    intervalIndex_ = labelMin;
    executionIndex_ = labelMin;
    offsetCount_ = 0;
    offsetPeriod_ = 0;

    switch (timeControl_)
    {
        case timeControls::timeStep:
            readIntervals<label>
            (
                intervalName,
                dict,
                nil(),
                true,
                1,
                intervalCounts_
            );
            break;

        case timeControls::writeTime:
        case timeControls::outputTime:
            readIntervals<label>
            (
                intervalName,
                dict,
                nil(),
                true,
                1,
                intervalCounts_
            );
            break;

        case timeControls::adjustableRunTime:
        case timeControls::runTime:
            readIntervals<scalar>
            (
                intervalName,
                dict,
                time_.userUnits(),
                intervalPeriods_
            );
            break;

        case timeControls::clockTime:
        case timeControls::cpuTime:
            readIntervals<scalar>
            (
                intervalName,
                dict,
                dimensions::time,
                intervalPeriods_
            );
            break;

        case timeControls::none:
            break;

        default:
            FatalErrorInFunction
                << "Undefined time control: "
                << timeControlNames_[timeControl_] << nl
                << exit(FatalError);
            break;
    }
}


bool Foam::timeControl::active() const
{
    return
        timeControl_ != timeControls::none
     && time_.value() >= startTime_ - 0.5*time_.deltaTValue()
     && time_.value() <= endTime_;
}


bool Foam::timeControl::execute()
{
    if (time_.writeTime()) writeTimeIndex_ ++;

    switch (timeControl_)
    {
        case timeControls::timeStep:
            return
                execute<label>
                (
                    intervalCounts_,
                    after(time_.value()),
                    0,
                    time_.timeIndex() - 1,
                    1,
                    offsetCount_
                );

        case timeControls::writeTime:
        case timeControls::outputTime:
            return
                execute<label>
                (
                    intervalCounts_,
                    after(time_.value()),
                    0,
                    writeTimeIndex_ - 1,
                    1,
                    offsetCount_
                );

        case timeControls::runTime:
        case timeControls::adjustableRunTime:
            return
                execute<scalar>
                (
                    intervalPeriods_,
                    after(time_.value()),
                    beginTime_,
                    time_.value() - time_.deltaTValue()/2,
                    time_.deltaTValue(),
                    offsetPeriod_
                );

        case timeControls::clockTime:
            return
                execute<scalar>
                (
                    intervalPeriods_,
                    after(time_.value()),
                    0,
                    returnReduce(label(time_.elapsedClockTime()), maxOp()),
                    0,
                    offsetPeriod_
                );

        case timeControls::cpuTime:
            return
                execute<scalar>
                (
                    intervalPeriods_,
                    after(time_.value()),
                    0,
                    returnReduce(label(time_.elapsedCpuTime()), maxOp()),
                    0,
                    offsetPeriod_
                );

        case timeControls::none:
            return false;

        default:
            FatalErrorInFunction
                << "Undefined time control: "
                << timeControlNames_[timeControl_] << nl
                << exit(FatalError);
            return false;
    }

    return false;
}


Foam::scalar Foam::timeControl::timeToNextAction() const
{
    if (time_.value() > endTime_) return vGreat;

    switch (timeControl_)
    {
        case timeControls::timeStep:
        case timeControls::writeTime:
        case timeControls::outputTime:
        case timeControls::runTime:
            return vGreat;

        case timeControls::adjustableRunTime:
        {
            const scalar tInterval =
                intervalIndex_ >= 0
              ? (
                    intervalIndex_ < intervalPeriods_.size() - 1
                  ? intervalTimes_[intervalIndex_ + 1]
                  + repeating(after(time_.value()))
                  : vGreat
                )
              : startTime_;

            scalar tNext =
                intervalIndex_ >= 0
              ? (
                    (executionIndex_ + 1)
                   *mag(intervalPeriods_[intervalIndex_])
                  + beginTime_ + offsetPeriod_
                )
              : vGreat;

            if (tNext > tInterval)
            {
                label intervalIndexStar = intervalIndex_;
                label executionIndexStar = executionIndex_;
                scalar offsetPeriodStar = offsetPeriod_;

                update<scalar>
                (
                    intervalPeriods_,
                    after(tInterval),
                    beginTime_,
                    tInterval - time_.deltaTValue()/2,
                    intervalIndexStar,
                    executionIndexStar,
                    offsetPeriodStar
                );

                tNext =
                    executionIndexStar < 0
                  ? tInterval
                  : (executionIndexStar + 1)
                   *mag(intervalPeriods_[intervalIndexStar])
                  + beginTime_ + offsetPeriodStar;
            }

            return tNext - time_.value();
        }

        case timeControls::clockTime:
        case timeControls::cpuTime:
        case timeControls::none:
            return vGreat;

        default:
            FatalErrorInFunction
                << "Undefined time control: "
                << timeControlNames_[timeControl_] << nl
                << exit(FatalError);
            return vGreat;
    }

    return vGreat;
}


// ************************************************************************* //
