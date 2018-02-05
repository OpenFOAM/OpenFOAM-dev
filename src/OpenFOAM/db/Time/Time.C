/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "Time.H"
#include "PstreamReduceOps.H"
#include "argList.H"
#include "IOdictionary.H"

#include <sstream>

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Time, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::Time::stopAtControls,
        4
    >::names[] =
    {
        "endTime",
        "noWriteNow",
        "writeNow",
        "nextWrite"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::Time::writeControls,
        5
    >::names[] =
    {
        "timeStep",
        "runTime",
        "adjustableRunTime",
        "clockTime",
        "cpuTime"
    };
}

const Foam::NamedEnum<Foam::Time::stopAtControls, 4>
    Foam::Time::stopAtControlNames_;

const Foam::NamedEnum<Foam::Time::writeControls, 5>
    Foam::Time::writeControlNames_;

Foam::Time::fmtflags Foam::Time::format_(Foam::Time::general);

int Foam::Time::precision_(6);

const int Foam::Time::maxPrecision_(3 - log10(small));

Foam::word Foam::Time::controlDictName("controlDict");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::Time::adjustDeltaT()
{
    scalar timeToNextWrite = max
    (
        0.0,
        (writeTimeIndex_ + 1)*writeInterval_ - (value() - startTime_)
    );

    timeToNextWrite = min(timeToNextWrite, functionObjects_.timeToNextWrite());

    scalar nSteps = timeToNextWrite/deltaT_;

    // For tiny deltaT the label can overflow!
    if (nSteps < labelMax)
    {
        label nStepsToNextWrite = label(nSteps + 0.5);

        scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

        // Control the increase of the time step to within a factor of 2
        // and the decrease within a factor of 5.
        if (newDeltaT >= deltaT_)
        {
            deltaT_ = min(newDeltaT, 2.0*deltaT_);
        }
        else
        {
            deltaT_ = max(newDeltaT, 0.2*deltaT_);
        }
    }
}


void Foam::Time::setControls()
{
    // default is to resume calculation from "latestTime"
    const word startFrom = controlDict_.lookupOrDefault<word>
    (
        "startFrom",
        "latestTime"
    );

    if (startFrom == "startTime")
    {
        controlDict_.lookup("startTime") >> startTime_;
    }
    else
    {
        // Search directory for valid time directories
        instantList timeDirs = findTimes(path(), constant());

        if (startFrom == "firstTime")
        {
            if (timeDirs.size())
            {
                if (timeDirs[0].name() == constant() && timeDirs.size() >= 2)
                {
                    startTime_ = timeDirs[1].value();
                }
                else
                {
                    startTime_ = timeDirs[0].value();
                }
            }
        }
        else if (startFrom == "latestTime")
        {
            if (timeDirs.size())
            {
                startTime_ = timeDirs.last().value();
            }
        }
        else
        {
            FatalIOErrorInFunction(controlDict_)
                << "expected startTime, firstTime or latestTime"
                << " found '" << startFrom << "'"
                << exit(FatalIOError);
        }
    }

    setTime(startTime_, 0);

    readDict();
    deltaTSave_ = deltaT_;
    deltaT0_ = deltaT_;

    // Check if time directory exists
    // If not increase time precision to see if it is formatted differently.
    if (!fileHandler().exists(timePath(), false))
    {
        int oldPrecision = precision_;
        int requiredPrecision = -1;
        bool found = false;
        word oldTime(timeName());
        for
        (
            precision_ = maxPrecision_;
            precision_ > oldPrecision;
            precision_--
        )
        {
            // Update the time formatting
            setTime(startTime_, 0);

            word newTime(timeName());
            if (newTime == oldTime)
            {
                break;
            }
            oldTime = newTime;

            // Check the existence of the time directory with the new format
            found = fileHandler().exists(timePath(), false);

            if (found)
            {
                requiredPrecision = precision_;
            }
        }

        if (requiredPrecision > 0)
        {
            // Update the time precision
            precision_ = requiredPrecision;

            // Update the time formatting
            setTime(startTime_, 0);

            WarningInFunction
                << "Increasing the timePrecision from " << oldPrecision
                << " to " << precision_
                << " to support the formatting of the current time directory "
                << timeName() << nl << endl;
        }
        else
        {
            // Could not find time directory so assume it is not present
            precision_ = oldPrecision;

            // Revert the time formatting
            setTime(startTime_, 0);
        }
    }

    if (Pstream::parRun())
    {
        scalar sumStartTime = startTime_;
        reduce(sumStartTime, sumOp<scalar>());
        if
        (
            mag(Pstream::nProcs()*startTime_ - sumStartTime)
          > Pstream::nProcs()*deltaT_/10.0
        )
        {
            FatalIOErrorInFunction(controlDict_)
                << "Start time is not the same for all processors" << nl
                << "processor " << Pstream::myProcNo() << " has startTime "
                << startTime_ << exit(FatalIOError);
        }
    }

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    // Read and set the deltaT only if time-step adjustment is active
    // otherwise use the deltaT from the controlDict
    if (controlDict_.lookupOrDefault<Switch>("adjustTimeStep", false))
    {
        if (timeDict.readIfPresent("deltaT", deltaT_))
        {
            deltaTSave_ = deltaT_;
            deltaT0_ = deltaT_;
        }
    }

    timeDict.readIfPresent("deltaT0", deltaT0_);

    if (timeDict.readIfPresent("index", startTimeIndex_))
    {
        timeIndex_ = startTimeIndex_;
    }


    // Check if values stored in time dictionary are consistent

    // 1. Based on time name
    bool checkValue = true;

    string storedTimeName;
    if (timeDict.readIfPresent("name", storedTimeName))
    {
        if (storedTimeName == timeName())
        {
            // Same time. No need to check stored value
            checkValue = false;
        }
    }

    // 2. Based on time value
    //    (consistent up to the current time writing precision so it won't
    //     trigger if we just change the write precision)
    if (checkValue)
    {
        scalar storedTimeValue;
        if (timeDict.readIfPresent("value", storedTimeValue))
        {
            word storedTimeName(timeName(storedTimeValue));

            if (storedTimeName != timeName())
            {
                IOWarningInFunction(timeDict)
                    << "Time read from time dictionary " << storedTimeName
                    << " differs from actual time " << timeName() << '.' << nl
                    << "    This may cause unexpected database behaviour."
                    << " If you are not interested" << nl
                    << "    in preserving time state delete"
                    << " the time dictionary."
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Time::Time
(
    const word& controlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(great),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_(*this, enableFunctionObjects)
{
    libs_.open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();

    // Time objects not registered so do like objectRegistry::checkIn ourselves.
    if (runTimeModifiable_)
    {
        // Monitor all files that controlDict depends on
        fileHandler().addWatches(controlDict_, controlDict_.files());
    }

    // Clear dependent files
    controlDict_.files().clear();
}


Foam::Time::Time
(
    const word& controlDictName,
    const argList& args,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        args.parRunControl().parRun(),
        args.rootPath(),
        args.globalCaseName(),
        args.caseName(),
        systemName,
        constantName
    ),

    objectRegistry(*this),

    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(great),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_
    (
        *this,
        argList::validOptions.found("withFunctionObjects")
      ? args.optionFound("withFunctionObjects")
      : !args.optionFound("noFunctionObjects")
    )
{
    libs_.open(controlDict_, "libs");

    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();

    // Time objects not registered so do like objectRegistry::checkIn ourselves.
    if (runTimeModifiable_)
    {
        // Monitor all files that controlDict depends on
        fileHandler().addWatches(controlDict_, controlDict_.files());
    }

    // Clear dependent files since not needed
    controlDict_.files().clear();
}


Foam::Time::Time
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dict
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(great),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),
    sigWriteNow_(true, *this),
    sigStopAtWriteNow_(true, *this),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_(*this, enableFunctionObjects)
{
    libs_.open(controlDict_, "libs");


    // Explicitly set read flags on objectRegistry so anything constructed
    // from it reads as well (e.g. fvSolution).
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    // Since could not construct regIOobject with setting:
    controlDict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;

    setControls();

    // Time objects not registered so do like objectRegistry::checkIn ourselves.
    if (runTimeModifiable_)
    {
        // Monitor all files that controlDict depends on
        fileHandler().addWatches(controlDict_, controlDict_.files());
    }

    // Clear dependent files since not needed
    controlDict_.files().clear();
}


Foam::Time::Time
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName,
    const bool enableFunctionObjects
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    libs_(),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(great),
    purgeWrite_(0),
    writeOnce_(false),
    subCycling_(false),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(false),

    functionObjects_(*this, enableFunctionObjects)
{
    libs_.open(controlDict_, "libs");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Time::~Time()
{
    forAllReverse(controlDict_.watchIndices(), i)
    {
        fileHandler().removeWatch(controlDict_.watchIndices()[i]);
    }

    // Destroy function objects first
    functionObjects_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::Time::timeName(const scalar t, const int precision)
{
    std::ostringstream buf;
    buf.setf(ios_base::fmtflags(format_), ios_base::floatfield);
    buf.precision(precision);
    buf << t;
    return buf.str();
}


Foam::word Foam::Time::timeName() const
{
    return dimensionedScalar::name();
}


Foam::instantList Foam::Time::times() const
{
    return findTimes(path(), constant());
}


Foam::word Foam::Time::findInstancePath
(
    const fileName& directory,
    const instant& t
) const
{
    // Simplified version: use findTimes (readDir + sort). The expensive
    // bit is the readDir, not the sorting. Tbd: avoid calling findInstancePath
    // from filePath.

    instantList timeDirs = findTimes(path(), constant());
    // Note:
    // - times will include constant (with value 0) as first element.
    //   For backwards compatibility make sure to find 0 in preference
    //   to constant.
    // - list is sorted so could use binary search

    forAllReverse(timeDirs, i)
    {
        if (t.equal(timeDirs[i].value()))
        {
            return timeDirs[i].name();
        }
    }

    return word::null;
}


Foam::word Foam::Time::findInstancePath(const instant& t) const
{
    return findInstancePath(path(), t);
}


Foam::instant Foam::Time::findClosestTime(const scalar t) const
{
    instantList timeDirs = findTimes(path(), constant());

    // There is only one time (likely "constant") so return it
    if (timeDirs.size() == 1)
    {
        return timeDirs[0];
    }

    if (t < timeDirs[1].value())
    {
        return timeDirs[1];
    }
    else if (t > timeDirs.last().value())
    {
        return timeDirs.last();
    }

    label nearestIndex = -1;
    scalar deltaT = great;

    for (label timei=1; timei < timeDirs.size(); ++timei)
    {
        scalar diff = mag(timeDirs[timei].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = timei;
        }
    }

    return timeDirs[nearestIndex];
}


Foam::label Foam::Time::findClosestTimeIndex
(
    const instantList& timeDirs,
    const scalar t,
    const word& constantName
)
{
    label nearestIndex = -1;
    scalar deltaT = great;

    forAll(timeDirs, timei)
    {
        if (timeDirs[timei].name() == constantName) continue;

        scalar diff = mag(timeDirs[timei].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = timei;
        }
    }

    return nearestIndex;
}


Foam::label Foam::Time::startTimeIndex() const
{
    return startTimeIndex_;
}


Foam::dimensionedScalar Foam::Time::startTime() const
{
    return dimensionedScalar("startTime", dimTime, startTime_);
}


Foam::dimensionedScalar Foam::Time::endTime() const
{
    return dimensionedScalar("endTime", dimTime, endTime_);
}


bool Foam::Time::running() const
{
    return value() < (endTime_ - 0.5*deltaT_);
}


bool Foam::Time::run() const
{
    bool running = this->running();

    if (!subCycling_)
    {
        if (!running && timeIndex_ != startTimeIndex_)
        {
            functionObjects_.execute();
            functionObjects_.end();
        }
    }

    if (running)
    {
        if (!subCycling_)
        {
            const_cast<Time&>(*this).readModifiedObjects();

            if (timeIndex_ == startTimeIndex_)
            {
                functionObjects_.start();
            }
            else
            {
                functionObjects_.execute();
            }
        }

        // Re-evaluate if running in case a function object has changed things
        running = this->running();
    }

    return running;
}


bool Foam::Time::loop()
{
    bool running = run();

    if (running)
    {
        operator++();
    }

    return running;
}


bool Foam::Time::end() const
{
    return value() > (endTime_ + 0.5*deltaT_);
}


bool Foam::Time::stopAt(const stopAtControls sa) const
{
    const bool changed = (stopAt_ != sa);
    stopAt_ = sa;

    // adjust endTime
    if (sa == saEndTime)
    {
        controlDict_.lookup("endTime") >> endTime_;
    }
    else
    {
        endTime_ = great;
    }
    return changed;
}


void Foam::Time::setTime(const Time& t)
{
    value() = t.value();
    dimensionedScalar::name() = t.dimensionedScalar::name();
    timeIndex_ = t.timeIndex_;
    fileHandler().setTime(*this);
}


void Foam::Time::setTime(const instant& inst, const label newIndex)
{
    value() = inst.value();
    dimensionedScalar::name() = inst.name();
    timeIndex_ = newIndex;

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    timeDict.readIfPresent("deltaT", deltaT_);
    timeDict.readIfPresent("deltaT0", deltaT0_);
    timeDict.readIfPresent("index", timeIndex_);
    fileHandler().setTime(*this);
}


void Foam::Time::setTime(const dimensionedScalar& newTime, const label newIndex)
{
    setTime(newTime.value(), newIndex);
}


void Foam::Time::setTime(const scalar newTime, const label newIndex)
{
    value() = newTime;
    dimensionedScalar::name() = timeName(timeToUserTime(newTime));
    timeIndex_ = newIndex;
    fileHandler().setTime(*this);
}


void Foam::Time::setEndTime(const dimensionedScalar& endTime)
{
    setEndTime(endTime.value());
}


void Foam::Time::setEndTime(const scalar endTime)
{
    endTime_ = endTime;
}


void Foam::Time::setDeltaT(const dimensionedScalar& deltaT)
{
    setDeltaT(deltaT.value());
}


void Foam::Time::setDeltaT(const scalar deltaT)
{
    setDeltaTNoAdjust(deltaT);

    functionObjects_.setTimeStep();

    if (writeControl_ == wcAdjustableRunTime)
    {
        adjustDeltaT();
    }
}


void Foam::Time::setDeltaTNoAdjust(const scalar deltaT)
{
    deltaT_ = deltaT;
    deltaTchanged_ = true;
}


Foam::TimeState Foam::Time::subCycle(const label nSubCycles)
{
    subCycling_ = true;
    prevTimeState_.set(new TimeState(*this));

    setTime(*this - deltaT(), (timeIndex() - 1)*nSubCycles);
    deltaT_ /= nSubCycles;
    deltaT0_ /= nSubCycles;
    deltaTSave_ = deltaT0_;

    return prevTimeState();
}


void Foam::Time::endSubCycle()
{
    if (subCycling_)
    {
        subCycling_ = false;
        TimeState::operator=(prevTimeState());
        prevTimeState_.clear();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Time& Foam::Time::operator+=(const dimensionedScalar& deltaT)
{
    return operator+=(deltaT.value());
}


Foam::Time& Foam::Time::operator+=(const scalar deltaT)
{
    setDeltaT(deltaT);
    return operator++();
}


Foam::Time& Foam::Time::operator++()
{
    deltaT0_ = deltaTSave_;
    deltaTSave_ = deltaT_;

    // Save old time value and name
    const scalar oldTimeValue = timeToUserTime(value());
    const word oldTimeName = dimensionedScalar::name();

    // Increment time
    setTime(value() + deltaT_, timeIndex_ + 1);

    if (!subCycling_)
    {
        // If the time is very close to zero reset to zero
        if (mag(value()) < 10*small*deltaT_)
        {
            setTime(0.0, timeIndex_);
        }

        if (sigStopAtWriteNow_.active() || sigWriteNow_.active())
        {
            // A signal might have been sent on one processor only
            // Reduce so all decide the same.

            label flag = 0;
            if (sigStopAtWriteNow_.active() && stopAt_ == saWriteNow)
            {
                flag += 1;
            }
            if (sigWriteNow_.active() && writeOnce_)
            {
                flag += 2;
            }
            reduce(flag, maxOp<label>());

            if (flag & 1)
            {
                stopAt_ = saWriteNow;
            }
            if (flag & 2)
            {
                writeOnce_ = true;
            }
        }

        writeTime_ = false;

        switch (writeControl_)
        {
            case wcTimeStep:
                writeTime_ = !(timeIndex_ % label(writeInterval_));
            break;

            case wcRunTime:
            case wcAdjustableRunTime:
            {
                label writeIndex = label
                (
                    ((value() - startTime_) + 0.5*deltaT_)
                  / writeInterval_
                );

                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcCpuTime:
            {
                label writeIndex = label
                (
                    returnReduce(elapsedCpuTime(), maxOp<double>())
                  / writeInterval_
                );
                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;

            case wcClockTime:
            {
                label writeIndex = label
                (
                    returnReduce(label(elapsedClockTime()), maxOp<label>())
                  / writeInterval_
                );
                if (writeIndex > writeTimeIndex_)
                {
                    writeTime_ = true;
                    writeTimeIndex_ = writeIndex;
                }
            }
            break;
        }


        // Check if endTime needs adjustment to stop at the next run()/end()
        if (!end())
        {
            if (stopAt_ == saNoWriteNow)
            {
                endTime_ = value();
            }
            else if (stopAt_ == saWriteNow)
            {
                endTime_ = value();
                writeTime_ = true;
            }
            else if (stopAt_ == saNextWrite && writeTime_ == true)
            {
                endTime_ = value();
            }
        }

        // Override writeTime if one-shot writing
        if (writeOnce_)
        {
            writeTime_ = true;
            writeOnce_ = false;
        }

        // Adjust the precision of the time directory name if necessary
        if (writeTime_)
        {
            // Tolerance used when testing time equivalence
            const scalar timeTol =
                max(min(pow(10.0, -precision_), 0.1*deltaT_), small);

            // User-time equivalent of deltaT
            const scalar userDeltaT = timeToUserTime(deltaT_);

            // Time value obtained by reading timeName
            scalar timeNameValue = -vGreat;

            // Check that new time representation differs from old one
            // reinterpretation of the word
            if
            (
                readScalar(dimensionedScalar::name().c_str(), timeNameValue)
             && (mag(timeNameValue - oldTimeValue - userDeltaT) > timeTol)
            )
            {
                int oldPrecision = precision_;
                while
                (
                    precision_ < maxPrecision_
                 && readScalar(dimensionedScalar::name().c_str(), timeNameValue)
                 && (mag(timeNameValue - oldTimeValue - userDeltaT) > timeTol)
                )
                {
                    precision_++;
                    setTime(value(), timeIndex());
                }

                if (precision_ != oldPrecision)
                {
                    WarningInFunction
                        << "Increased the timePrecision from " << oldPrecision
                        << " to " << precision_
                        << " to distinguish between timeNames at time "
                        << dimensionedScalar::name()
                        << endl;

                    if (precision_ == maxPrecision_)
                    {
                        // Reached maxPrecision limit
                        WarningInFunction
                            << "Current time name " << dimensionedScalar::name()
                            << nl
                            << "    The maximum time precision has been reached"
                               " which might result in overwriting previous"
                               " results."
                            << endl;
                    }

                    // Check if round-off error caused time-reversal
                    scalar oldTimeNameValue = -vGreat;
                    if
                    (
                        readScalar(oldTimeName.c_str(), oldTimeNameValue)
                     && (
                            sign(timeNameValue - oldTimeNameValue)
                         != sign(deltaT_)
                        )
                    )
                    {
                        WarningInFunction
                            << "Current time name " << dimensionedScalar::name()
                            << " is set to an instance prior to the "
                               "previous one "
                            << oldTimeName << nl
                            << "    This might result in temporal "
                               "discontinuities."
                            << endl;
                    }
                }
            }
        }
    }

    return *this;
}


Foam::Time& Foam::Time::operator++(int)
{
    return operator++();
}


// ************************************************************************* //
