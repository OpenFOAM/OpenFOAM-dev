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
#include "Pstream.H"
#include "simpleObjectRegistry.H"
#include "dimensionedConstants.H"
#include "IOdictionary.H"
#include "fileOperation.H"

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
        InfoHeader
            << "Overriding DebugSwitches according to " << controlDict_.name()
            << endl;

        simpleObjectRegistry& objects = debug::debugObjects();
        const dictionary& localSettings = controlDict_.subDict("DebugSwitches");
        forAllConstIter(dictionary, localSettings, iter)
        {
            const word& name = iter().keyword();

            simpleObjectRegistryEntry* objPtr = objects.lookupPtr(name);

            if (objPtr)
            {
                InfoHeader << "    " << iter() << endl;

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
        InfoHeader
            << "Overriding OptimisationSwitches according to "
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
                InfoHeader << "    " << iter() << endl;

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


        // Handle fileHandler override explicitly since interacts with
        // local dictionary monitoring.
        word fileHandlerName;
        if
        (
            localSettings.readIfPresent("fileHandler", fileHandlerName)
         && fileHandler().type() != fileHandlerName
        )
        {
            // Remove the old watches since destroying the file
            fileNameList oldWatchedFiles(controlDict_.watchIndices());
            forAllReverse(controlDict_.watchIndices(), i)
            {
                label watchi = controlDict_.watchIndices()[i];
                oldWatchedFiles[i] = fileHandler().getFile(watchi);
                fileHandler().removeWatch(watchi);
            }
            controlDict_.watchIndices().clear();

            // Installing the new handler
            InfoHeader
                << "Overriding fileHandler to " << fileHandlerName << endl;

            autoPtr<fileOperation> handler
            (
                fileOperation::New
                (
                    fileHandlerName,
                    true
                )
            );
            Foam::fileHandler(handler);

            // Reinstall old watches
            fileHandler().addWatches(controlDict_, oldWatchedFiles);
        }
    }


    // DimensionedConstants. Handled as a special case since both e.g.
    // the 'unitSet' might be changed and the individual values
    if (controlDict_.found("DimensionedConstants"))
    {
        InfoHeader
            << "Overriding DimensionedConstants according to "
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

                if (writeInfoHeader)
                {

                    Info<< "    ";
                    objects[i]->writeData(Info);
                    Info<< endl;
                }
            }
        }
    }


    // Dimension sets
    if (controlDict_.found("DimensionSets"))
    {
        InfoHeader
            << "Overriding DimensionSets according to "
            << controlDict_.name() << endl;

        dictionary dict(Foam::dimensionSystems());
        dict.merge(controlDict_.subDict("DimensionSets"));

        simpleObjectRegistry& objects = debug::dimensionSetObjects();

        simpleObjectRegistryEntry* objPtr = objects.lookupPtr("DimensionSets");

        if (objPtr)
        {
            InfoHeader << controlDict_.subDict("DimensionSets") << endl;

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

    if (controlDict_.readIfPresent("writeInterval", writeInterval_))
    {
        if (writeControl_ == wcTimeStep && label(writeInterval_) < 1)
        {
            FatalIOErrorInFunction(controlDict_)
                << "writeInterval < 1 for writeControl timeStep"
                << exit(FatalIOError);
        }
    }
    else
    {
        controlDict_.lookup("writeFrequency") >> writeInterval_;
    }


    if (oldWriteInterval != writeInterval_)
    {
        switch (writeControl_)
        {
            case wcRunTime:
            case wcAdjustableRunTime:
                // Recalculate writeTimeIndex_ to be in units of current
                // writeInterval.
                writeTimeIndex_ = label
                (
                    writeTimeIndex_
                  * oldWriteInterval
                  / writeInterval_
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
            WarningInFunction
                << "invalid value for purgeWrite " << purgeWrite_
                << ", should be >= 0, setting to 0"
                << endl;

            purgeWrite_ = 0;
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
            WarningInFunction
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
            endTime_ = great;
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

        if
        (
            writeFormat_ == IOstream::BINARY
         && writeCompression_ == IOstream::COMPRESSED
        )
        {
            IOWarningInFunction(controlDict_)
                << "Selecting compressed binary is inefficient and ineffective"
                   ", resetting to uncompressed binary"
                << endl;

            writeCompression_ = IOstream::UNCOMPRESSED;
        }
    }

    controlDict_.readIfPresent("graphFormat", graphFormat_);
    controlDict_.readIfPresent("runTimeModifiable", runTimeModifiable_);




    if (!runTimeModifiable_ && controlDict_.watchIndices().size())
    {
        forAllReverse(controlDict_.watchIndices(), i)
        {
            fileHandler().removeWatch(controlDict_.watchIndices()[i]);
        }
        controlDict_.watchIndices().clear();
    }
}


bool Foam::Time::read()
{
    if (controlDict_.regIOobject::read())
    {
        readDict();

        if (runTimeModifiable_)
        {
            // For IOdictionary the call to regIOobject::read() would have
            // already updated all the watchIndices via the addWatch but
            // controlDict_ is an unwatchedIOdictionary so will only have
            // stored the dependencies as files.
            fileHandler().addWatches(controlDict_, controlDict_.files());
        }
        controlDict_.files().clear();

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
        fileHandler().updateStates
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

            if (runTimeModifiable_)
            {
                // For IOdictionary the call to regIOobject::read() would have
                // already updated all the watchIndices via the addWatch but
                // controlDict_ is an unwatchedIOdictionary so will only have
                // stored the dependencies as files.

                fileHandler().addWatches(controlDict_, controlDict_.files());
            }
            controlDict_.files().clear();
        }

        bool registryModified = objectRegistry::modified();

        if (registryModified)
        {
            objectRegistry::readModifiedObjects();
        }
    }
}


bool Foam::Time::writeTimeDict() const
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

    timeDict.add("value", timeName(timeToUserTime(value()), maxPrecision_));
    timeDict.add("name", string(tmName));
    timeDict.add("index", timeIndex_);
    timeDict.add("deltaT", timeToUserTime(deltaT_));
    timeDict.add("deltaT0", timeToUserTime(deltaT0_));

    return timeDict.regIOobject::writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED,
        true
    );
}


bool Foam::Time::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    if (writeTime())
    {
        bool writeOK = writeTimeDict();

        if (writeOK)
        {
            writeOK = objectRegistry::writeObject(fmt, ver, cmp, valid);
        }

        if (writeOK)
        {
            // Does the writeTime trigger purging?
            if (writeTime_ && purgeWrite_)
            {
                if
                (
                    previousWriteTimes_.size() == 0
                 || previousWriteTimes_.top() != timeName()
                )
                {
                    previousWriteTimes_.push(timeName());
                }

                while (previousWriteTimes_.size() > purgeWrite_)
                {
                    fileHandler().rmDir
                    (
                        fileHandler().filePath
                        (
                            objectRegistry::path(previousWriteTimes_.pop())
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
    writeTime_ = true;
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
