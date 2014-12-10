/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "Pstream.H"
#include "simpleObjectRegistry.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Time::readDict()
{
    word application;
    if (controlDict_.readIfPresent("application", application))
    {
        // Do not override if already set so external application can override
        setEnv("FOAM_APPLICATION", application, false);
    }


    // Check for local switches and settings
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Debug switches
    if (controlDict_.found("DebugSwitches"))
    {
        Info<< "Overriding DebugSwitches according to " << controlDict_.name()
            << endl;

        simpleObjectRegistry& objects = debug::debugObjects();
        const dictionary& localSettings = controlDict_.subDict("DebugSwitches");
        forAllConstIter(dictionary, localSettings, iter)
        {
            const word& name = iter().keyword();

            simpleObjectRegistryEntry* objPtr = objects.lookupPtr(name);

            if (objPtr)
            {
                Info<< "    " << iter() << endl;

                const List<simpleRegIOobject*>& objects = *objPtr;

                if (iter().isDict())
                {
                    forAll(objects, i)
                    {
                        OStringStream os(IOstream::ASCII);
                        os  << iter().dict();
                        IStringStream is(os.str());
                        objects[i]->readData(is);
                    }
                }
                else
                {
                    forAll(objects, i)
                    {
                        objects[i]->readData(iter().stream());
                    }
                }
            }
        }
    }

    // Optimisation Switches
    if (controlDict_.found("OptimisationSwitches"))
    {
        Info<< "Overriding OptimisationSwitches according to "
            << controlDict_.name() << endl;

        simpleObjectRegistry& objects = debug::optimisationObjects();
        const dictionary& localSettings = controlDict_.subDict
        (
            "OptimisationSwitches"
        );
        forAllConstIter(dictionary, localSettings, iter)
        {
            const word& name = iter().keyword();

            simpleObjectRegistryEntry* objPtr = objects.lookupPtr(name);

            if (objPtr)
            {
                Info<< "    " << iter() << endl;

                const List<simpleRegIOobject*>& objects = *objPtr;

                if (iter().isDict())
                {
                    forAll(objects, i)
                    {
                        OStringStream os(IOstream::ASCII);
                        os  << iter().dict();
                        IStringStream is(os.str());
                        objects[i]->readData(is);
                    }
                }
                else
                {
                    forAll(objects, i)
                    {
                        objects[i]->readData(iter().stream());
                    }
                }
            }
        }
    }


    // DimensionedConstants. Handled as a special case since both e.g.
    // the 'unitSet' might be changed and the individual values
    if (controlDict_.found("DimensionedConstants"))
    {
        Info<< "Overriding DimensionedConstants according to "
            << controlDict_.name() << endl;

        // Change in-memory
        dimensionedConstants().merge
        (
            controlDict_.subDict("DimensionedConstants")
        );


        simpleObjectRegistry& objects = debug::dimensionedConstantObjects();

        IStringStream dummyIs("");

        forAllConstIter(simpleObjectRegistry, objects, iter)
        {
            const List<simpleRegIOobject*>& objects = *iter;

            forAll(objects, i)
            {
                objects[i]->readData(dummyIs);

                Info<< "    ";
                objects[i]->writeData(Info);
                Info<< endl;
            }
        }
    }


    // Dimension sets
    if (controlDict_.found("DimensionSets"))
    {
        Info<< "Overriding DimensionSets according to "
            << controlDict_.name() << endl;

        dictionary dict(Foam::dimensionSystems());
        dict.merge(controlDict_.subDict("DimensionSets"));

        simpleObjectRegistry& objects = debug::dimensionSetObjects();

        simpleObjectRegistryEntry* objPtr = objects.lookupPtr("DimensionSets");

        if (objPtr)
        {
            Info<< controlDict_.subDict("DimensionSets") << endl;

            const List<simpleRegIOobject*>& objects = *objPtr;

            forAll(objects, i)
            {
                OStringStream os(IOstream::ASCII);
                os  << dict;
                IStringStream is(os.str());
                objects[i]->readData(is);
            }
        }
    }


    if (!deltaTchanged_)
    {
        deltaT_ = readScalar(controlDict_.lookup("deltaT"));
    }

    if (controlDict_.found("writeControl"))
    {
        writeControl_ = writeControlNames_.read
        (
            controlDict_.lookup("writeControl")
        );
    }

    scalar oldWriteInterval = writeInterval_;
    scalar oldSecondaryWriteInterval = secondaryWriteInterval_;

    if (controlDict_.readIfPresent("writeInterval", writeInterval_))
    {
        if (writeControl_ == wcTimeStep && label(writeInterval_) < 1)
        {
            FatalIOErrorIn("Time::readDict()", controlDict_)
                << "writeInterval < 1 for writeControl timeStep"
                << exit(FatalIOError);
        }
    }
    else
    {
        controlDict_.lookup("writeFrequency") >> writeInterval_;
    }


    // Additional writing
    if (controlDict_.found("secondaryWriteControl"))
    {
        secondaryWriteControl_ = writeControlNames_.read
        (
            controlDict_.lookup("secondaryWriteControl")
        );

        if
        (
            controlDict_.readIfPresent
            (
                "secondaryWriteInterval",
                secondaryWriteInterval_
            )
        )
        {
            if
            (
                secondaryWriteControl_ == wcTimeStep
             && label(secondaryWriteInterval_) < 1
            )
            {
                FatalIOErrorIn("Time::readDict()", controlDict_)
                    << "secondaryWriteInterval < 1"
                    << " for secondaryWriteControl timeStep"
                    << exit(FatalIOError);
            }
        }
        else
        {
            controlDict_.lookup("secondaryWriteFrequency")
                >> secondaryWriteInterval_;
        }
    }



    if (oldWriteInterval != writeInterval_)
    {
        switch (writeControl_)
        {
            case wcRunTime:
            case wcAdjustableRunTime:
                // Recalculate outputTimeIndex_ to be in units of current
                // writeInterval.
                outputTimeIndex_ = label
                (
                    outputTimeIndex_
                  * oldWriteInterval
                  / writeInterval_
                );
            break;

            default:
            break;
        }
    }
    if (oldSecondaryWriteInterval != secondaryWriteInterval_)
    {
        switch (secondaryWriteControl_)
        {
            case wcRunTime:
            case wcAdjustableRunTime:
                // Recalculate secondaryOutputTimeIndex_ to be in units of
                // current writeInterval.
                secondaryOutputTimeIndex_ = label
                (
                    secondaryOutputTimeIndex_
                  * oldSecondaryWriteInterval
                  / secondaryWriteInterval_
                );
            break;

            default:
            break;
        }
    }

    if (controlDict_.readIfPresent("purgeWrite", purgeWrite_))
    {
        if (purgeWrite_ < 0)
        {
            WarningIn("Time::readDict()")
                << "invalid value for purgeWrite " << purgeWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            purgeWrite_ = 0;
        }
    }

    if (controlDict_.readIfPresent("secondaryPurgeWrite", secondaryPurgeWrite_))
    {
        if (secondaryPurgeWrite_ < 0)
        {
            WarningIn("Time::readDict()")
                << "invalid value for secondaryPurgeWrite "
                << secondaryPurgeWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            secondaryPurgeWrite_ = 0;
        }
    }

    if (controlDict_.found("timeFormat"))
    {
        const word formatName(controlDict_.lookup("timeFormat"));

        if (formatName == "general")
        {
            format_ = general;
        }
        else if (formatName == "fixed")
        {
            format_ = fixed;
        }
        else if (formatName == "scientific")
        {
            format_ = scientific;
        }
        else
        {
            WarningIn("Time::readDict()")
                << "unsupported time format " << formatName
                << endl;
        }
    }

    controlDict_.readIfPresent("timePrecision", precision_);

    // stopAt at 'endTime' or a specified value
    // if nothing is specified, the endTime is zero
    if (controlDict_.found("stopAt"))
    {
        stopAt_ = stopAtControlNames_.read(controlDict_.lookup("stopAt"));

        if (stopAt_ == saEndTime)
        {
            controlDict_.lookup("endTime") >> endTime_;
        }
        else
        {
            endTime_ = GREAT;
        }
    }
    else if (!controlDict_.readIfPresent("endTime", endTime_))
    {
        endTime_ = 0;
    }

    dimensionedScalar::name() = timeName(value());

    if (controlDict_.found("writeVersion"))
    {
        writeVersion_ = IOstream::versionNumber
        (
            controlDict_.lookup("writeVersion")
        );
    }

    if (controlDict_.found("writeFormat"))
    {
        writeFormat_ = IOstream::formatEnum
        (
            controlDict_.lookup("writeFormat")
        );
    }

    if (controlDict_.found("writePrecision"))
    {
        IOstream::defaultPrecision
        (
            readUint(controlDict_.lookup("writePrecision"))
        );

        Sout.precision(IOstream::defaultPrecision());
        Serr.precision(IOstream::defaultPrecision());

        Pout.precision(IOstream::defaultPrecision());
        Perr.precision(IOstream::defaultPrecision());

        FatalError().precision(IOstream::defaultPrecision());
        FatalIOError.error::operator()().precision
        (
            IOstream::defaultPrecision()
        );
    }

    if (controlDict_.found("writeCompression"))
    {
        writeCompression_ = IOstream::compressionEnum
        (
            controlDict_.lookup("writeCompression")
        );
    }

    controlDict_.readIfPresent("graphFormat", graphFormat_);
    controlDict_.readIfPresent("runTimeModifiable", runTimeModifiable_);

    if (!runTimeModifiable_ && controlDict_.watchIndex() != -1)
    {
        removeWatch(controlDict_.watchIndex());
        controlDict_.watchIndex() = -1;
    }
}


bool Foam::Time::read()
{
    if (controlDict_.regIOobject::read())
    {
        readDict();
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::Time::readModifiedObjects()
{
    if (runTimeModifiable_)
    {
        // Get state of all monitored objects (=registered objects with a
        // valid filePath).
        // Note: requires same ordering in objectRegistries on different
        // processors!
        monitorPtr_().updateStates
        (
            (
                regIOobject::fileModificationChecking == inotifyMaster
             || regIOobject::fileModificationChecking == timeStampMaster
            ),
            Pstream::parRun()
        );

        // Time handling is special since controlDict_ is the one dictionary
        // that is not registered to any database.

        if (controlDict_.readIfModified())
        {
            readDict();
            functionObjects_.read();
        }

        bool registryModified = objectRegistry::modified();

        if (registryModified)
        {
            objectRegistry::readModifiedObjects();
        }
    }
}


bool Foam::Time::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    if (outputTime())
    {
        const word tmName(timeName());

        IOdictionary timeDict
        (
            IOobject
            (
                "time",
                tmName,
                "uniform",
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        timeDict.add("value", timeToUserTime(value()));
        timeDict.add("name", string(tmName));
        timeDict.add("index", timeIndex_);
        timeDict.add("deltaT", timeToUserTime(deltaT_));
        timeDict.add("deltaT0", timeToUserTime(deltaT0_));

        timeDict.regIOobject::writeObject(fmt, ver, cmp);
        bool writeOK = objectRegistry::writeObject(fmt, ver, cmp);

        if (writeOK)
        {
            // Does primary or secondary time trigger purging?
            // Note that primary times can only be purged by primary
            // purging. Secondary times can be purged by either primary
            // or secondary purging.
            if (primaryOutputTime_ && purgeWrite_)
            {
                previousOutputTimes_.push(tmName);

                while (previousOutputTimes_.size() > purgeWrite_)
                {
                    rmDir(objectRegistry::path(previousOutputTimes_.pop()));
                }
            }
            if
            (
               !primaryOutputTime_
             && secondaryOutputTime_
             && secondaryPurgeWrite_
            )
            {
                // Writing due to secondary
                previousSecondaryOutputTimes_.push(tmName);

                while
                (
                    previousSecondaryOutputTimes_.size()
                  > secondaryPurgeWrite_
                )
                {
                    rmDir
                    (
                        objectRegistry::path
                        (
                            previousSecondaryOutputTimes_.pop()
                        )
                    );
                }
            }
        }

        return writeOK;
    }
    else
    {
        return false;
    }
}


bool Foam::Time::writeNow()
{
    primaryOutputTime_ = true;
    outputTime_ = true;
    return write();
}


bool Foam::Time::writeAndEnd()
{
    stopAt_  = saWriteNow;
    endTime_ = value();

    return writeNow();
}


void Foam::Time::writeOnce()
{
    writeOnce_ = true;
}


// ************************************************************************* //
