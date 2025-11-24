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

Description
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "etcFiles.H"
#include "Ostream.H"
#include "demandDrivenData.H"
#include "IOobject.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

dictionary* configDictPtr_(nullptr);

dictionary* debugSwitchesPtr_(nullptr);
dictionary* infoSwitchesPtr_(nullptr);
dictionary* optimisationSwitchesPtr_(nullptr);

dictionary* debugDefaultSwitchesPtr_(nullptr);
dictionary* infoDefaultSwitchesPtr_(nullptr);
dictionary* optimisationDefaultSwitchesPtr_(nullptr);

dictionary& debugDefaultSwitches()
{
    if (!debugDefaultSwitchesPtr_)
    {
        debugDefaultSwitchesPtr_ = new dictionary();
    }

    return *debugDefaultSwitchesPtr_;
}

dictionary& infoDefaultSwitches()
{
    if (!infoDefaultSwitchesPtr_)
    {
        infoDefaultSwitchesPtr_ = new dictionary();
    }

    return *infoDefaultSwitchesPtr_;
}

dictionary& optimisationDefaultSwitches()
{
    if (!optimisationDefaultSwitchesPtr_)
    {
        optimisationDefaultSwitchesPtr_ = new dictionary();
    }

    return *optimisationDefaultSwitchesPtr_;
}


// To ensure configDictPtr_ is deleted at the end of the run
class deleteControlDictPtr
{
public:

    deleteControlDictPtr()
    {}

    ~deleteControlDictPtr()
    {
        deleteDemandDrivenData(debugDefaultSwitchesPtr_);
        deleteDemandDrivenData(infoDefaultSwitchesPtr_);
        deleteDemandDrivenData(optimisationDefaultSwitchesPtr_);

        debugSwitchesPtr_ = nullptr;
        infoSwitchesPtr_ = nullptr;
        optimisationSwitchesPtr_ = nullptr;

        deleteDemandDrivenData(configDictPtr_);
    }
};

deleteControlDictPtr deleteControlDictPtr_;
//! \endcond


} // End namespace debug
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::debug::configDict()
{
    if (!configDictPtr_)
    {
        string configDictString(getEnv("FOAM_CONFIGDICT"));
        if (!configDictString.empty())
        {
            // Read from environment
            IStringStream is(configDictString);
            configDictPtr_ = new dictionary(is);
        }
        else
        {
            fileNameList configDictFiles = findEtcFiles("configDict", true);
            configDictPtr_ = new dictionary();
            forAllReverse(configDictFiles, cdfi)
            {
                IFstream ifs(configDictFiles[cdfi]);

                if (!ifs.good())
                {
                    SafeFatalIOErrorInFunction
                    (
                        ifs,
                        "Cannot open configDict"
                    );
                }
                configDictPtr_->merge(dictionary(ifs));
            }

            // Check for legacy etc controlDict files
            // for backwards compatibility
            fileNameList controlDictFiles = findEtcFiles("controlDict", false);

            if (controlDictFiles.size())
            {
                cout<< "--> FOAM Warning: legacy controlDict"
                       " configuration files found:" << std::endl;

                forAll(controlDictFiles, i)
                {
                    cout<< "    " << controlDictFiles[i] << std::endl;
                }

                cout<< "    Please rename these files controlDict -> configDict"
                    << std::endl;
            }

            forAllReverse(controlDictFiles, cdfi)
            {
                IFstream ifs(controlDictFiles[cdfi]);

                if (!ifs.good())
                {
                    SafeFatalIOErrorInFunction
                    (
                        ifs,
                        "Cannot open controlDict"
                    );
                }
                configDictPtr_->merge(dictionary(ifs));
            }


        }

        IFstream ifs("system/configDict");
        if (ifs.good())
        {
            entry::disableFunctionEntries = true;
            configDictPtr_->merge(dictionary(ifs));
        }
    }

    return *configDictPtr_;
}


Foam::dictionary& Foam::debug::switchSet
(
    const char* subDictName,
    dictionary*& subDictPtr
)
{
    if (!subDictPtr)
    {
        entry* ePtr = configDict().lookupEntryPtr
        (
            subDictName, false, false
        );

        if (!ePtr || !ePtr->isDict())
        {
            cerr<< "debug::switchSet(const char*, dictionary*&):\n"
                << "    Cannot find " <<  subDictName << " in dictionary "
                << configDict().name().c_str()
                << std::endl << std::endl;

            ::exit(1);
        }

        subDictPtr = &ePtr->dict();
    }

    return *subDictPtr;
}


Foam::dictionary& Foam::debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


Foam::dictionary& Foam::debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


Foam::dictionary& Foam::debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


int Foam::debug::debugSwitch(const char* name, const int defaultValue)
{
    if
    (
        debugDefaultSwitches().found(name)
     && debugDefaultSwitches().lookup<int>(name) != defaultValue
    )
    {
        FatalErrorInFunction
            << "Multiple defaults set for debug switch " << name
            << exit(FatalError);
    }

    debugDefaultSwitches().set(name, defaultValue);

    return debugSwitches().lookupOrAddDefault(name, defaultValue);
}


int Foam::debug::infoSwitch(const char* name, const int defaultValue)
{
    if
    (
        infoDefaultSwitches().found(name)
     && infoDefaultSwitches().lookup<int>(name) != defaultValue
    )
    {
        FatalErrorInFunction
            << "Multiple defaults set for info switch " << name
            << exit(FatalError);
    }

    infoDefaultSwitches().set(name, defaultValue);

    return infoSwitches().lookupOrAddDefault(name, defaultValue);
}


int Foam::debug::optimisationSwitch(const char* name, const int defaultValue)
{
    if
    (
        optimisationDefaultSwitches().found(name)
     && optimisationDefaultSwitches().lookup<int>(name) != defaultValue
    )
    {
        FatalErrorInFunction
            << "Multiple defaults set for optimisation switch " << name
            << exit(FatalError);
    }

    optimisationDefaultSwitches().set(name, defaultValue);

    return optimisationSwitches().lookupOrAddDefault(name, defaultValue);
}


float Foam::debug::floatOptimisationSwitch
(
    const char* name,
    const float defaultValue
)
{
    if
    (
        optimisationDefaultSwitches().found(name)
     && optimisationDefaultSwitches().lookup<float>(name) != defaultValue
    )
    {
        FatalErrorInFunction
            << exit(FatalError);
    }

    optimisationDefaultSwitches().set(name, defaultValue);

    return optimisationSwitches().lookupOrAddDefault(name, defaultValue);
}


const Foam::word Foam::debug::wordOptimisationSwitch
(
    const char* name,
    const word& defaultValue
)
{
    if
    (
        optimisationDefaultSwitches().found(name)
     && optimisationDefaultSwitches().lookup<word>(name) != defaultValue
    )
    {
        FatalErrorInFunction
            << exit(FatalError);
    }

    optimisationDefaultSwitches().set(name, defaultValue);

    return optimisationSwitches().lookupOrAddDefault(name, defaultValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void listSwitches
(
    const wordList& debugSwitches,
    const wordList& infoSwitches,
    const wordList& optSwitches,
    const bool unset
)
{
    if (unset)
    {
        fileNameList configDictFiles = findEtcFiles("configDict", true);
        dictionary configDict;
        forAllReverse(configDictFiles, cdfi)
        {
            configDict.merge(dictionary(IFstream(configDictFiles[cdfi])()));
        }

        wordHashSet configDictDebug
        (
            configDict.subDict("DebugSwitches").sortedToc()
        );

        wordHashSet configDictInfo
        (
            configDict.subDict("InfoSwitches").sortedToc()
        );

        wordHashSet configDictOpt
        (
            configDict.subDict("OptimisationSwitches").sortedToc()
        );


        IOobject::writeDivider(Info);

        wordHashSet hashset;
        hashset = debugSwitches;
        hashset -= configDictDebug;
        Info<< "Unset DebugSwitches" << hashset.sortedToc() << endl;

        hashset = infoSwitches;
        hashset -= configDictInfo;
        Info<< "Unset InfoSwitches" << hashset.sortedToc() << endl;

        hashset = optSwitches;
        hashset -= configDictOpt;
        Info<< "Unset OptimisationSwitches" << hashset.sortedToc() << endl;
    }
    else
    {
        IOobject::writeDivider(Info);
        Info<< "DebugSwitches" << debugSwitches << endl;
        Info<< "InfoSwitches" << infoSwitches << endl;
        Info<< "OptimisationSwitches" << optSwitches << endl;
    }
}


void listSwitches
(
    const word& name,
    const dictionary& switches,
    const dictionary& defaultSwitches
)
{
    wordHashSet defaultSet;
    wordHashSet nonDefaultSet;
    wordHashSet noDefaultSet;

    forAllConstIter(dictionary, switches, iter)
    {
        const word& name = iter().keyword();

        const bool hasDefault = defaultSwitches.found(name);

        const bool isDefault =
            hasDefault
         && defaultSwitches.lookupEntry(name, false, false) == iter();

        if (hasDefault)
        {
            if (isDefault)
            {
                defaultSet.insert(name);
            }
            else
            {
                nonDefaultSet.insert(name);
            }
        }
        else
        {
            noDefaultSet.insert(name);
        }
    }

    auto print = [&](const char* heading, const wordList& names)
    {
        Info<< indent << "// " << heading << endl;

        forAll(names, i)
        {
            Info<< switches.lookupEntry(names[i], false, false);
        }
    };

    Info<< name << endl
        << token::BEGIN_BLOCK << endl << incrIndent;

    print
    (
        "Switches with default values",
        defaultSet.sortedToc()
    );
    Info<< nl;
    print
    (
        "Switches with non-default values",
        nonDefaultSet.sortedToc()
    );
    Info<< nl;
    print
    (
        "Switches without defaults",
        noDefaultSet.sortedToc()
    );

    Info<< decrIndent << token::END_BLOCK << endl;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::debug::listSwitches()
{
    listSwitches
    (
        "DebugSwitches",
        debug::debugSwitches(),
        debug::debugDefaultSwitches()
    );
    Info<< endl;

    listSwitches
    (
        "InfoSwitches",
        debug::infoSwitches(),
        debug::infoDefaultSwitches()
    );
    Info<< endl;

    listSwitches
    (
        "OptimisationSwitches",
        debug::optimisationSwitches(),
        debug::optimisationDefaultSwitches()
    );
}


// ************************************************************************* //
