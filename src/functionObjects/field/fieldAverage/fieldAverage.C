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

#include "fieldAverage.H"
#include "fieldAverageItem.H"
#include "timeIOdictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldAverage, 0);
    addToRunTimeSelectionTable(functionObject, fieldAverage, dictionary);
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldAverage::baseType,
    2
>::names[] = { "iteration", "time"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldAverage::baseType,
    2
> Foam::functionObjects::fieldAverage::baseTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldAverage::initialise()
{
    // Initialise any unset times
    forAll(totalTime_, fieldi)
    {
        if (totalTime_[fieldi] < 0)
        {
            totalTime_[fieldi] = 0;
        }
    }

    // Initialise mean fields
    forAll(faItems_, fieldi)
    {
        initialiseMeanField<scalar>(fieldi);
        initialiseMeanField<vector>(fieldi);
        initialiseMeanField<sphericalTensor>(fieldi);
        initialiseMeanField<symmTensor>(fieldi);
        initialiseMeanField<tensor>(fieldi);
    }
    forAll(faItems_, fieldi)
    {
        initialisePrime2MeanField<scalar, scalar>(fieldi);
        initialisePrime2MeanField<vector, symmTensor>(fieldi);
    }
}


void Foam::functionObjects::fieldAverage::restart()
{
    Log << "    Restarting averaging at time " << time_.name() << nl;

    // Clear the times
    totalIter_ = 0;
    totalTime_ = -1;

    // Clear mean fields
    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            if (obr_.found(faItems_[i].meanFieldName()))
            {
                obr_.checkOut(*obr_[faItems_[i].meanFieldName()]);
            }
        }

        if (faItems_[i].prime2Mean())
        {
            if (obr_.found(faItems_[i].prime2MeanFieldName()))
            {
                obr_.checkOut(*obr_[faItems_[i].prime2MeanFieldName()]);
            }
        }
    }

    // Re-create any mean fields
    initialise();

    // Ensure first averaging works unconditionally
    prevTimeIndex_ = -1;
}


void Foam::functionObjects::fieldAverage::calcAverages()
{
    Log << type() << " " << name() << ":" << nl;

    const label currentTimeIndex = time_.timeIndex();
    const scalar currentTime = time_.value();

    if (prevTimeIndex_ == currentTimeIndex)
    {
        return;
    }

    prevTimeIndex_ = currentTimeIndex;

    if (periodicRestart_ && currentTime > restartPeriod_*periodIndex_)
    {
        restart();
        periodIndex_++;
    }
    else
    {
        initialise();
    }

    // Increment the time and iteration totals
    forAll(faItems_, fieldi)
    {
        totalIter_[fieldi]++;
        totalTime_[fieldi] += time_.deltaTValue();
    }

    Log << "    Calculating averages" << nl;

    addMeanSqrToPrime2Mean<scalar, scalar>();
    addMeanSqrToPrime2Mean<vector, symmTensor>();

    calculateMeanFields<scalar>();
    calculateMeanFields<vector>();
    calculateMeanFields<sphericalTensor>();
    calculateMeanFields<symmTensor>();
    calculateMeanFields<tensor>();

    calculatePrime2MeanFields<scalar, scalar>();
    calculatePrime2MeanFields<vector, symmTensor>();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::writeAverages() const
{
    Log << type() << " " << name() << ":" << nl
        << "    Writing average fields" << endl;

    writeFields<scalar>();
    writeFields<vector>();
    writeFields<sphericalTensor>();
    writeFields<symmTensor>();
    writeFields<tensor>();

    timeIOdictionary propsDict
    (
        IOobject
        (
            name() + "Properties",
            time_.name(),
            "uniform",
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(faItems_, fieldi)
    {
        const word& fieldName = faItems_[fieldi].fieldName();
        propsDict.add(fieldName, dictionary());
        propsDict.subDict(fieldName).add("totalIter", totalIter_[fieldi]);
        propsDict.subDict(fieldName).add("totalTime", totalTime_[fieldi]);
    }

    propsDict.regIOobject::write();

    Log << endl;
}


void Foam::functionObjects::fieldAverage::read
(
    const dictionary& dict,
    const bool construct
)
{
    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    dict.readIfPresent("restartOnOutput", restartOnOutput_);
    dict.readIfPresent("periodicRestart", periodicRestart_);

    if (periodicRestart_)
    {
        dict.lookup("restartPeriod") >> restartPeriod_;
    }

    mean_ = dict.lookupOrDefault<Switch>("mean", true);
    prime2Mean_ = dict.lookupOrDefault<Switch>("prime2Mean", false);
    base_ = baseTypeNames_[dict.lookupOrDefault<word>("base", "time")];
    window_ = dict.lookupOrDefault<scalar>("window", -1);
    windowName_ = dict.lookupOrDefault<word>("windowName", "");

    if (construct)
    {
        // First read of a run. Look for properties dict and read total
        // iter/time values from it if available.

        faItems_.clear();
        totalIter_.clear();
        totalTime_.clear();

        faItems_ =
            PtrList<fieldAverageItem>
            (
                dict.lookup("fields"),
                fieldAverageItem::iNew(*this)
            );
        totalIter_.setSize(faItems_.size(), 0);
        totalTime_.setSize(faItems_.size(), -1);

        typeIOobject<timeIOdictionary> propsDictHeader
        (
            name() + "Properties",
            time_.name(),
            "uniform",
            obr_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        );

        dictionary propsDict;
        if
        (
            !restartOnRestart_
         && !restartOnOutput_
         && propsDictHeader.headerOk()
        )
        {
            propsDict = timeIOdictionary(propsDictHeader);
        }

        bool first = true;
        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            if (!propsDict.found(fieldName))
            {
                if (first)
                {
                    Log << "    Starting averaging for fields:" << nl;
                    first = false;
                }
                Log << "        " << fieldName << nl;
            }
        }

        first = true;
        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            if (propsDict.found(fieldName))
            {
                if (first)
                {
                    Log << "    Restarting averaging for fields:" << nl;
                    first = false;
                }
                const dictionary& fieldPropsDict = propsDict.subDict(fieldName);
                totalIter_[fieldi] = fieldPropsDict.lookup<label>("totalIter");
                totalTime_[fieldi] = fieldPropsDict.lookup<scalar>("totalTime");
                Log << "        " << fieldName
                    << " iters = " << totalIter_[fieldi]
                    << " time = " << totalTime_[fieldi] << nl;
            }
        }
    }
    else
    {
        // Re-read during a run. Read the items and copy the per-field total
        // iter/times from the old items to the new.

        PtrList<fieldAverageItem> faItems0;
        List<label> totalIter0;
        List<scalar> totalTime0;
        faItems0.transfer(faItems_);
        totalIter0.transfer(totalIter_);
        totalTime0.transfer(totalTime_);

        faItems_ =
            PtrList<fieldAverageItem>
            (
                dict.lookup("fields"),
                fieldAverageItem::iNew(*this)
            );
        totalIter_.resize(faItems_.size(), 0);
        totalTime_.resize(faItems_.size(), -1);

        // Map from field to old-field index
        labelList fieldiFieldi0s(faItems_.size(), -1);
        forAll(faItems_, fieldi)
        {
            const word& fieldName = faItems_[fieldi].fieldName();
            forAll(faItems0, fieldi0)
            {
                if (faItems0[fieldi0].fieldName() == fieldName)
                {
                    fieldiFieldi0s[fieldi] = fieldi0;
                    break;
                }
            }
        }

        bool first = true;
        forAll(faItems_, fieldi)
        {
            if (fieldiFieldi0s[fieldi] == -1)
            {
                if (first)
                {
                    Log << "    Starting averaging for fields:" << nl;
                    first = true;
                }
                Log << "        " << faItems_[fieldi].fieldName() << nl;
            }
        }

        first = true;
        forAll(faItems_, fieldi)
        {
            if (fieldiFieldi0s[fieldi] != -1)
            {
                if (first)
                {
                    Log << "    Continuing averaging for fields:" << nl;
                    first = true;
                }
                totalIter_[fieldi] = totalIter0[fieldiFieldi0s[fieldi]];
                totalTime_[fieldi] = totalTime0[fieldiFieldi0s[fieldi]];
                Log << "        " << faItems_[fieldi].fieldName()
                    << " iters = " << totalIter_[fieldi]
                    << " time = " << totalTime_[fieldi] << nl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::fieldAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    prevTimeIndex_(-1),
    restartOnRestart_(false),
    restartOnOutput_(false),
    periodicRestart_(false),
    restartPeriod_(great),
    base_(baseType::iter),
    window_(-1.0),
    windowName_(""),
    faItems_(),
    totalIter_(),
    totalTime_(),
    periodIndex_(1)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name << ":" << nl;

    read(dict, true);

    // Read any available mean fields
    forAll(faItems_, fieldi)
    {
        readMeanField<scalar>(fieldi);
        readMeanField<vector>(fieldi);
        readMeanField<sphericalTensor>(fieldi);
        readMeanField<symmTensor>(fieldi);
        readMeanField<tensor>(fieldi);
    }
    forAll(faItems_, fieldi)
    {
        readPrime2MeanField<scalar, scalar>(fieldi);
        readPrime2MeanField<vector, symmTensor>(fieldi);
    }

    Log << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverage::~fieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldAverage::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    read(dict, false);

    Log << endl;

    return true;
}


Foam::wordList Foam::functionObjects::fieldAverage::fields() const
{
    wordList fields(faItems_.size());

    forAll(faItems_, fieldi)
    {
        fields[fieldi] = faItems_[fieldi].fieldName();
    }

    return fields;
}


bool Foam::functionObjects::fieldAverage::execute()
{
    if (functionObject::postProcess)
    {
        WarningInFunction
            << "fieldAverage is not supported with the foamPostProcess utility"
            << endl;

        return false;
    }

    calcAverages();

    return true;
}


bool Foam::functionObjects::fieldAverage::write()
{
    writeAverages();

    if (restartOnOutput_)
    {
        restart();
    }

    return true;
}


// ************************************************************************* //
