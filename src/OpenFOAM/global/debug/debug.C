/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

//! \cond ignoreDocumentation
//- Skip documentation : local scope only

dictionary* controlDictPtr_(nullptr);

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


// To ensure controlDictPtr_ is deleted at the end of the run
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

        deleteDemandDrivenData(controlDictPtr_);
    }
};

deleteControlDictPtr deleteControlDictPtr_;
//! \endcond


} // End namespace debug
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::debug::controlDict()
{
    if (!controlDictPtr_)
    {
        string controlDictString(getEnv("FOAM_CONTROLDICT"));
        if (!controlDictString.empty())
        {
            // Read from environment
            IStringStream is(controlDictString);
            controlDictPtr_ = new dictionary(is);
        }
        else
        {
            fileNameList controlDictFiles = findEtcFiles("controlDict", true);
            controlDictPtr_ = new dictionary();
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
                controlDictPtr_->merge(dictionary(ifs));
            }
        }

        IFstream ifs("system/controlDict");
        if (ifs.good())
        {
            entry::disableFunctionEntries = true;
            controlDictPtr_->merge(dictionary(ifs));
        }
    }

    return *controlDictPtr_;
}


Foam::dictionary& Foam::debug::switchSet
(
    const char* subDictName,
    dictionary*& subDictPtr
)
{
    if (!subDictPtr)
    {
        entry* ePtr = controlDict().lookupEntryPtr
        (
            subDictName, false, false
        );

        if (!ePtr || !ePtr->isDict())
        {
            cerr<< "debug::switchSet(const char*, dictionary*&):\n"
                << "    Cannot find " <<  subDictName << " in dictionary "
                << controlDict().name().c_str()
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

    return debugSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
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

    return infoSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
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

    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
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

    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
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

    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
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
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        dictionary controlDict;
        forAllReverse(controlDictFiles, cdfi)
        {
            controlDict.merge(dictionary(IFstream(controlDictFiles[cdfi])()));
        }

        wordHashSet controlDictDebug
        (
            controlDict.subDict("DebugSwitches").sortedToc()
        );

        wordHashSet controlDictInfo
        (
            controlDict.subDict("InfoSwitches").sortedToc()
        );

        wordHashSet controlDictOpt
        (
            controlDict.subDict("OptimisationSwitches").sortedToc()
        );


        IOobject::writeDivider(Info);

        wordHashSet hashset;
        hashset = debugSwitches;
        hashset -= controlDictDebug;
        Info<< "Unset DebugSwitches" << hashset.sortedToc() << endl;

        hashset = infoSwitches;
        hashset -= controlDictInfo;
        Info<< "Unset InfoSwitches" << hashset.sortedToc() << endl;

        hashset = optSwitches;
        hashset -= controlDictOpt;
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
