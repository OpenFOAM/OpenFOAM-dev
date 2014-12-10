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

Description
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "Ostream.H"
#include "demandDrivenData.H"
#include "simpleObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

//! \cond ignoreDocumentation - local scope
dictionary* controlDictPtr_(NULL);
dictionary* debugSwitchesPtr_(NULL);
dictionary* infoSwitchesPtr_(NULL);
dictionary* optimisationSwitchesPtr_(NULL);

// Debug switch read and write callback tables.
simpleObjectRegistry* debugObjectsPtr_(NULL);
simpleObjectRegistry* infoObjectsPtr_(NULL);
simpleObjectRegistry* optimisationObjectsPtr_(NULL);
simpleObjectRegistry* dimensionSetObjectsPtr_(NULL);
simpleObjectRegistry* dimensionedConstantObjectsPtr_(NULL);


// to ensure controlDictPtr_ is deleted at the end of the run
class deleteControlDictPtr
{
public:

    deleteControlDictPtr()
    {}

    ~deleteControlDictPtr()
    {
        deleteDemandDrivenData(debugObjectsPtr_);
        deleteDemandDrivenData(infoObjectsPtr_);
        deleteDemandDrivenData(optimisationObjectsPtr_);
        deleteDemandDrivenData(dimensionSetObjectsPtr_);
        deleteDemandDrivenData(dimensionedConstantObjectsPtr_);

        debugSwitchesPtr_ = NULL;
        infoSwitchesPtr_ = NULL;
        optimisationSwitchesPtr_ = NULL;
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
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        controlDictPtr_ = new dictionary();
        forAllReverse(controlDictFiles, cdfi)
        {
            IFstream ifs(controlDictFiles[cdfi]);

            if (!ifs.good())
            {
                SafeFatalIOErrorIn
                (
                    "debug::controlDict()",
                    ifs,
                    "Cannot open controlDict"
                );
            }
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
    return debugSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::infoSwitch(const char* name, const int defaultValue)
{
    return infoSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::optimisationSwitch(const char* name, const int defaultValue)
{
    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


void Foam::debug::addDebugObject(const char* name, simpleRegIOobject* obj)
{
    simpleObjectRegistryEntry* ptr = debugObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        debugObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addInfoObject(const char* name, simpleRegIOobject* obj)
{
    simpleObjectRegistryEntry* ptr = infoObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        infoObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addOptimisationObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = optimisationObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        optimisationObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addDimensionSetObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = dimensionSetObjects().lookupPtr(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        dimensionSetObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


void Foam::debug::addDimensionedConstantObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = dimensionedConstantObjects().lookupPtr
    (
        name
    );
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        dimensionedConstantObjects().append
        (
            name,
            new simpleObjectRegistryEntry
            (
                List<simpleRegIOobject*>(1, obj)
            )
        );
    }
}


Foam::simpleObjectRegistry& Foam::debug::debugObjects()
{
    if (!debugObjectsPtr_)
    {
        debugObjectsPtr_ = new simpleObjectRegistry(1000);
    }

    return *debugObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::infoObjects()
{
    if (!infoObjectsPtr_)
    {
        infoObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *infoObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::optimisationObjects()
{
    if (!optimisationObjectsPtr_)
    {
        optimisationObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *optimisationObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::dimensionSetObjects()
{
    if (!dimensionSetObjectsPtr_)
    {
        dimensionSetObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *dimensionSetObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::dimensionedConstantObjects()
{
    if (!dimensionedConstantObjectsPtr_)
    {
        dimensionedConstantObjectsPtr_ = new simpleObjectRegistry(100);
    }

    return *dimensionedConstantObjectsPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
