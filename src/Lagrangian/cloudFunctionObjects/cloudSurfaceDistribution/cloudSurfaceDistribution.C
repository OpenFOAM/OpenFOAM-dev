/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloudSurfaceDistribution.H"
#include "cloud.H"
#include "faceSet.H"
#include "functionObject.H"
#include "hashedWordList.H"
#include "OSspecific.H"
#include "timeIOdictionary.H"
#include "unintegrable.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudSurfaceDistribution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cloudSurfaceDistribution,
        dictionary
    );
}
}

const Foam::NamedEnum
<
    Foam::functionObjects::cloudSurfaceDistribution::selectionType,
    3
>
Foam::functionObjects::cloudSurfaceDistribution::selectionTypeNames
{"faceZone", "faceSet", "patch"};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::cloudSurfaceDistribution::readFields
(
    const dictionary& dict,
    const word& key,
    const wordList& defaultValue
)
{
    const bool haveFields = dict.found(key + "s");
    const bool haveField = dict.found(key);

    if (notNull(defaultValue) && !haveFields && !haveField)
    {
        return defaultValue;
    }

    if (haveFields == haveField)
    {
        FatalIOErrorInFunction(dict)
            << "keywords " << key << "s and " << key << " both "
            << (haveFields ? "" : "un") << "defined in "
            << "dictionary " << dict.name()
            << exit(FatalIOError);
    }

    if (haveFields)
    {
        return dict.lookup<wordList>(key + "s");
    }
    else
    {
        return wordList({dict.lookup<word>(key)});
    }
}


Foam::functionObjects::cloudSurfaceDistribution::selectionType
Foam::functionObjects::cloudSurfaceDistribution::readSelectionType
(
    const dictionary& dict
)
{
    // Use the "select" entry if it is there
    if (dict.found("select"))
    {
        return selectionTypeNames.read(dict.lookup("select"));
    }

    // If there is only one of the possible keys (i.e., "faceZone", "faceSet",
    // and "patch") present in the dictionary then infer the selection type
    // from that key
    const NamedEnum<selectionType, 3>::namesType& keys = selectionTypeNames;
    label keyi = -1;
    forAll(keys, i)
    {
         if (!dict.found(keys[i])) continue;

         keyi = keyi == -1 ? i : -labelMax;
    }
    if (keyi >= 0)
    {
        return selectionTypeNames[keys[keyi]];
    }

    // Raise an error requesting a "select" entry
    return selectionTypeNames.read(dict.lookup("select"));
}


const char* const*
Foam::functionObjects::cloudSurfaceDistribution::componentNames
(
    const label fieldi
) const
{
    #define FIELD_RANGE_COMPONENT_NAMES(Type, GeoField)                        \
        if (mesh().foundObject<GeoField<Type>>(fields_[fieldi]))               \
        {                                                                      \
            return pTraits<Type>::componentNames;                              \
        }

    FOR_ALL_FIELD_TYPES(FIELD_RANGE_COMPONENT_NAMES, LagrangianField)
    FOR_ALL_FIELD_TYPES(FIELD_RANGE_COMPONENT_NAMES, LagrangianDynamicField)
    FOR_ALL_FIELD_TYPES(FIELD_RANGE_COMPONENT_NAMES, LagrangianInternalField)

    return nullptr;
}


void Foam::functionObjects::cloudSurfaceDistribution::readCoeffs
(
    const dictionary& dict,
    const bool props
)
{
    // Re-read the fields
    const hashedWordList oldFields(fields_);
    fields_ = readFields(dict, "field");

    // Determine if any fields are being continued
    bool continued = false;
    forAll(fields_, fieldi)
    {
        continued = continued || oldFields.found(fields_[fieldi]);
    }

    // If we are continuing then the weight fields must not have changed
    const wordList oldWeightFields(weightFields_);
    weightFields_ = readFields(dict, "weightField", wordList());
    if (continued && weightFields_ != oldWeightFields)
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change weight fields at run-time"
            << exit(FatalIOError);
    }

    // If we are continuing then the surface must not have changed
    const selectionType oldSelectionType = selectionType_;
    selectionType_ = readSelectionType(dict);
    const word oldSelectionName = selectionName_;
    selectionName_ = dict.lookup<word>(selectionTypeNames[selectionType_]);
    if
    (
        selectionType_ != oldSelectionType
     || selectionName_ != oldSelectionName
    )
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the selected surface at run-time"
            << exit(FatalError);
    }

    // If we are continuing then the number of bins must not have changed
    const label oldNBins = nBins_;
    nBins_ = dict.lookup<label>("nBins");
    if (continued && nBins_ != oldNBins)
    {
        FatalIOErrorInFunction(dict)
            << "Cannot change the number of bins at run-time"
            << exit(FatalIOError);
    }

    // Re-read the formatter. This can change whenever.
    if (!props)
    {
        formatter_ = setWriter::New(dict.lookup("setFormat"), dict);
    }

    // If we have a new list of fields then we need to move any continuing
    // data to its new location
    if (fields_ != oldFields)
    {
        List<List<scalarField>> oldSums;
        oldSums.transfer(sums_);
        sums_.resize(fields_.size());

        List<List<Pair<scalar>>> oldRanges;
        oldRanges.transfer(ranges_);
        ranges_.resize(fields_.size());

        labelList oldNSamples;
        oldNSamples.transfer(nSamples_);
        nSamples_.resize(fields_.size());

        forAll(fields_, fieldi)
        {
            if (!oldFields.found(fields_[fieldi])) continue;

            const label oldFieldi = oldFields[fields_[fieldi]];

            sums_[fieldi].transfer(oldSums[oldFieldi]);
            ranges_[fieldi].transfer(oldRanges[oldFieldi]);
            nSamples_[fieldi] = oldNSamples[oldFieldi];
        }
    }
}


Foam::IOobject Foam::functionObjects::cloudSurfaceDistribution::propsDictIo
(
    const IOobject::readOption r
) const
{
    return
        IOobject
        (
            name() + "Properties",
            mesh().time().name(),
            "uniform",
            mesh(),
            r,
            IOobject::NO_WRITE,
            false
        );
}


Foam::boolList Foam::functionObjects::cloudSurfaceDistribution::selected
(
    const LagrangianSubScalarSubField& fraction
) const
{
    const Foam::cloud& c = cloud();

    const polyBoundaryMesh& pbm = c.mesh().mesh().boundaryMesh();

    const label patchi =
        max
        (
            static_cast<label>(fraction.mesh().group())
          - static_cast<label>(LagrangianGroup::onPatchZero),
            -1
        );

    const List<LagrangianState> states =
        c.mesh().changing()
      ? List<LagrangianState>(fraction.mesh().sub(c.mesh().states()))
      : patchi != -1
      ? List<LagrangianState>(fraction.size(), LagrangianState::onPatchZero)
      : List<LagrangianState>(fraction.size(), LagrangianState::none);

    const SubField<label> facei = fraction.mesh().sub(c.mesh().facei());

    boolList result(facei.size(), false);

    switch (selectionType_)
    {
        case selectionType::faceZone:
        {
            const faceZone& z = c.mesh().mesh().faceZones()[selectionName_];
            forAll(facei, subi)
            {
                if (states[subi] <= LagrangianState::inCell) continue;
                result = z.lookupMap().found(facei[subi]);
            }
            break;
        }
        case selectionType::faceSet:
        {
            forAll(facei, subi)
            {
                if (states[subi] <= LagrangianState::inCell) continue;
                result = selectionSet_.found(facei[subi]);
            }
            break;
        }
        case selectionType::patch:
        {
            result = patchi == pbm[selectionName_].index();
            break;
        }
    }

    return result;
}


template<template<class> class GeoField>
bool Foam::functionObjects::cloudSurfaceDistribution::multiplyWeight
(
    const LagrangianSubMesh& subMesh,
    const label weightFieldi,
    scalarField& weight
) const
{
    if (!mesh().foundObject<GeoField<scalar>>(weightFields_[weightFieldi]))
        return false;

    const GeoField<scalar>& w =
        mesh().lookupObject<GeoField<scalar>>(weightFields_[weightFieldi]);

    weight *= subMesh.sub(w.primitiveField());

    return true;
}


template<template<class> class GeoField, class Type>
bool Foam::functionObjects::cloudSurfaceDistribution::addField
(
    const LagrangianSubMesh& subMesh,
    const boolList& selected,
    const scalarField& weight,
    const label fieldi
)
{
    if (!mesh().foundObject<GeoField<Type>>(fields_[fieldi])) return false;

    const GeoField<Type>& field =
        mesh().lookupObject<GeoField<Type>>(fields_[fieldi]);

    // Allocate the data if needed
    if (sums_[fieldi].empty())
    {
        sums_[fieldi].resize(pTraits<Type>::nComponents);
        ranges_[fieldi].resize(pTraits<Type>::nComponents);

        for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
        {
            sums_[fieldi][d].resize(nBins_ + 1, scalar(0));
            ranges_[fieldi][d] = Pair<scalar>(vGreat, -vGreat);
        }

        nSamples_[fieldi] = 0;
    }

    // Consider each component in turn
    for (direction d = 0; d < pTraits<Type>::nComponents; ++ d)
    {
        Pair<scalar>& range = ranges_[fieldi][d];
        const Pair<scalar> range0 = ranges_[fieldi][d];

        // Expand the range to include the new elements if necessary
        bool rangeHasChanged = false;

        forAll(weight, subi)
        {
            if (!selected[subi]) continue;

            const scalar x = component(field[subMesh.start() + subi], d);

            if (x < range.first())
            {
                range.first() = x;
                rangeHasChanged = true;
            }

            if (x > range.second())
            {
                range.second() = x;
                rangeHasChanged = true;
            }
        }

        // Synchronise
        reduce(rangeHasChanged, orOp<bool>());
        if (rangeHasChanged)
        {
            reduce(range.first(), minOp<scalar>());
            reduce(range.second(), maxOp<scalar>());
        }

        // If the range was expanded then re-sample the old sums onto the new
        // discretisation of property space
        if (rangeHasChanged && nSamples_[fieldi])
        {
            scalarList sum0;
            sum0.transfer(sums_[fieldi][d]);
            sums_[fieldi][d].resize(nBins_ + 1, scalar(0));

            forAll(sum0, nodei0)
            {
                const scalar x =
                    (1 - scalar(nodei0)/nBins_)*range0.first()
                  + scalar(nodei0)/nBins_*range0.second();
                const scalar f =
                    (x - range.first())
                   /max(range.second() - range.first(), rootVSmall);
                const label bini = min(max(floor(f*nBins_), 0), nBins_ - 1);
                const scalar g = f*nBins_ - scalar(bini);

                sums_[fieldi][d][bini] += sum0[nodei0]*(1 - g);
                sums_[fieldi][d][bini + 1] += sum0[nodei0]*g;
            }
        }

        // Add the new elements to the bins
        forAll(weight, subi)
        {
            if (!selected[subi]) continue;

            const scalar x = component(field[subMesh.start() + subi], d);
            const scalar f =
                (x - range.first())
               /max(range.second() - range.first(), rootVSmall);
            const label bini = min(max(floor(f*nBins_), 0), nBins_ - 1);
            const scalar g = f*nBins_ - scalar(bini);

            sums_[fieldi][d][bini] += weight[subi]*(1 - g);
            sums_[fieldi][d][bini + 1] += weight[subi]*g;
        }
    }

    // Update the number of samples
    forAll(weight, subi)
    {
        if (!selected[subi]) continue;

        nSamples_[fieldi] ++;
    }

    return true;
}


void Foam::functionObjects::cloudSurfaceDistribution::writeDistribution
(
    const word& fieldName,
    const word& componentName,
    const scalarField& x,
    const scalarField& PDF,
    const scalarField& CDF
) const
{
    if (!Pstream::master()) return;

    const fileName& outputPath =
        time_.globalPath()
       /writeFile::outputPrefix
       /(
            mesh().mesh().name() != polyMesh::defaultRegion
          ? mesh().mesh().name()
          : word::null
        )
       /name()
       /time_.name();

    mkDir(outputPath);

    const word fieldComponentName =
        fieldName
      + (componentName.empty() ? "" : "_")
      + componentName;

    formatter_->write
    (
        outputPath,
        fieldComponentName,
        coordSet(true, fieldComponentName, x),
        "PDF", PDF,
        "CDF", CDF
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudSurfaceDistribution::cloudSurfaceDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    LagrangianMeshFunctionObject(name, runTime, dict, cloud::typeName),
    cloudFunctionObject
    (
        static_cast<const LagrangianMeshFunctionObject&>(*this)
    ),
    fields_(),
    weightFields_(readFields(dict, "weightField", wordList())),
    selectionType_(readSelectionType(dict)),
    selectionName_(dict.lookup<word>(selectionTypeNames[selectionType_])),
    selectionSet_
    (
        selectionType_ == selectionType::faceSet
      ? faceSet(mesh().mesh(), selectionName_)
      : labelHashSet()
    ),
    nBins_(dict.lookup<label>("nBins")),
    formatter_(),
    sums_(),
    ranges_(),
    nSamples_()
{
    // Restore from saved properties dictionary
    typeIOobject<timeIOdictionary> propsDictIo
    (
        this->propsDictIo(IOobject::MUST_READ)
    );
    if (propsDictIo.headerOk())
    {
        timeIOdictionary propsDict(propsDictIo);
        readCoeffs(propsDict, true);
        sums_ = propsDict.lookup<List<List<scalarField>>>("sums");
        ranges_ = propsDict.lookup<List<List<Pair<scalar>>>>("ranges");
        nSamples_ = propsDict.lookup<labelList>("nSamples");
    }

    // Read the coefficients from the function object dictionary
    readCoeffs(dict, false);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudSurfaceDistribution::~cloudSurfaceDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cloudSurfaceDistribution::read
(
    const dictionary& dict
)
{
    if (LagrangianMeshFunctionObject::read(dict))
    {
        readCoeffs(dict, false);
        return true;
    }
    else
    {
        return false;
    }
}


Foam::wordList Foam::functionObjects::cloudSurfaceDistribution::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::cloudSurfaceDistribution::executeAtStart() const
{
    return false;
}


bool Foam::functionObjects::cloudSurfaceDistribution::execute()
{
    return true;
}


void Foam::functionObjects::cloudSurfaceDistribution::preCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{}


void Foam::functionObjects::cloudSurfaceDistribution::postCrossFaces
(
    const LagrangianSubScalarSubField& fraction
)
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    if
    (
        selectionType_ == selectionType::patch
     && subMesh.group() < LagrangianGroup::onPatchZero
    ) return;

    // Identify which particles are to be included
    const boolList selected(this->selected(fraction));

    // Construct the weights
    scalarField weight(subMesh.size(), scalar(1));
    forAll(weightFields_, weightFieldi)
    {
        #define MULTIPLY_WEIGHT(GeoField) \
            && !multiplyWeight<GeoField>(subMesh, weightFieldi, weight)

        if
        (
            true
            MULTIPLY_WEIGHT(LagrangianField)
            MULTIPLY_WEIGHT(LagrangianDynamicField)
            MULTIPLY_WEIGHT(LagrangianInternalField)
        )
        {
            FatalErrorInFunction
                << "Weight field " << weightFields_[weightFieldi]
                << " was not found" << exit(FatalError);
        }
    }

    // Add to the sums
    forAll(fields_, fieldi)
    {
        #define ADD_FIELD(Type, GeoField) \
            && !addField<GeoField, Type>(subMesh, selected, weight, fieldi)

        if
        (
            true
            FOR_ALL_FIELD_TYPES(ADD_FIELD, LagrangianField)
            FOR_ALL_FIELD_TYPES(ADD_FIELD, LagrangianDynamicField)
            FOR_ALL_FIELD_TYPES(ADD_FIELD, LagrangianInternalField)
        )
        {
            cannotFindObject(fields_[fieldi]);
        }
    }
}


bool Foam::functionObjects::cloudSurfaceDistribution::write()
{
    forAll(sums_, fieldi)
    {
        // Quit if this PDF has nothing in it yet
        if (returnReduce(nSamples_[fieldi], sumOp<label>()) == 0) continue;

        // Get the component names for this field
        const char* const* componentNames = this->componentNames(fieldi);

        // Write each component's PDF
        forAll(sums_[fieldi], d)
        {
            const scalarList& sum = sums_[fieldi][d];
            const Pair<scalar>& range = ranges_[fieldi][d];

            // Write a single point if the distribution is uniform
            if (range.first() == range.second())
            {
                writeDistribution
                (
                    fields_[fieldi],
                    componentNames[d],
                    scalarField(1, range.first()),
                    scalarField(1, vGreat),
                    scalarField(1, scalar(1))
                );
                continue;
            }

            // Construct the limits of the bins
            scalarField x(nBins_ + 1);
            forAll(sum, nodei)
            {
                const scalar f = scalar(nodei)/nBins_;
                x[nodei] = (1 - f)*range.first() + f*range.second();
            }

            // Convert the sum to a PDF by synchronising, normalising and
            // correcting the ends (which have half the sample space of the
            // interior points)
            scalarField PDF(sum);
            Pstream::listCombineGather(PDF, plusEqOp<scalar>());
            Pstream::listCombineScatter(PDF);
            PDF /= Foam::sum(PDF)*(range.second() - range.first())/nBins_;
            PDF.first() *= 2;
            PDF.last() *= 2;

            // Write
            writeDistribution
            (
                fields_[fieldi],
                componentNames[d],
                x,
                PDF,
                distributions::unintegrable::integrate(x, PDF)()
            );
        }
    }

    if (!mesh().time().writeTime()) return true;

    timeIOdictionary propsDict(propsDictIo(IOobject::NO_READ));
    propsDict.add("select", selectionTypeNames[selectionType_]);
    propsDict.add(selectionTypeNames[selectionType_], selectionName_);
    propsDict.add("fields", fields_);
    propsDict.add("nBins", nBins_);
    propsDict.add("sums", sums_);
    propsDict.add("ranges", ranges_);
    propsDict.add("nSamples", nSamples_);
    propsDict.regIOobject::write();

    return true;
}


bool Foam::functionObjects::cloudSurfaceDistribution::clear()
{
    return true;
}


void Foam::functionObjects::cloudSurfaceDistribution::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (selectionType_ == selectionType::faceSet)
    {
        selectionSet_ = faceSet(mesh().mesh(), selectionName_);
    }
}


void Foam::functionObjects::cloudSurfaceDistribution::mapMesh
(
    const polyMeshMap&
)
{
    if (selectionType_ == selectionType::faceSet)
    {
        selectionSet_ = faceSet(mesh().mesh(), selectionName_);
    }
}


void Foam::functionObjects::cloudSurfaceDistribution::distribute
(
    const polyDistributionMap&
)
{
    if (selectionType_ == selectionType::faceSet)
    {
        selectionSet_ = faceSet(mesh().mesh(), selectionName_);
    }
}


// ************************************************************************* //
