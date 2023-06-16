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

#include "Pstream.H"
#include "functionObjectList.H"
#include "streamlines.H"
#include "streamlinesCloud.H"
#include "ReadFields.H"
#include "meshSearch.H"
#include "sampledSet.H"
#include "globalIndex.H"
#include "distributionMap.H"
#include "interpolationCellPoint.H"
#include "PatchTools.H"
#include "writeFile.H"
#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void gatherAndFlatten(DynamicField<Type>& field)
{
    List<List<Type>> gatheredField(Pstream::nProcs());
    gatheredField[Pstream::myProcNo()] = field;
    Pstream::gatherList(gatheredField);

    field =
        ListListOps::combine<List<Type>>
        (
            gatheredField,
            accessOp<List<Type>>()
        );
}

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
        NamedEnum<functionObjects::streamlines::trackDirection, 3>::names[] =
        {"forward", "backward", "both"};

    namespace functionObjects
    {
        defineTypeNameAndDebug(streamlines, 0);
        addToRunTimeSelectionTable(functionObject, streamlines, dictionary);

        const NamedEnum<streamlines::trackDirection, 3>
            streamlines::trackDirectionNames_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::streamlines::streamlines
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),
    nSubCycle_(0)
{
    read(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::streamlines::~streamlines()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::streamlines::read(const dictionary& dict)
{
    if (dict != dict_)
    {
        dict_ = dict;
    }

    Info<< type() << " " << name() << ":" << nl;

    dict.lookup("fields") >> fields_;

    UName_ = dict.lookupOrDefault("U", word("U"));

    writeAge_ = dict.lookupOrDefault<Switch>("writeAge", true);

    trackDirection_ = trackDirectionNames_[word(dict.lookup("direction"))];

    trackOutside_ = dict.lookupOrDefault<Switch>("outside", false);

    dict.lookup("lifeTime") >> lifeTime_;

    if (lifeTime_ < 1)
    {
        FatalErrorInFunction
            << "Illegal value " << lifeTime_ << " for lifeTime"
            << exit(FatalError);
    }

    bool subCycling = dict.found("nSubCycle");
    bool fixedLength = dict.found("trackLength");
    if (subCycling && fixedLength)
    {
        FatalIOErrorInFunction(dict)
            << "Cannot both specify automatic time stepping (through '"
            << "nSubCycle' specification) and fixed track length (through '"
            << "trackLength')"
            << exit(FatalIOError);
    }
    if (subCycling)
    {
        nSubCycle_ = max(dict.lookup<scalar>("nSubCycle"), 1);
        trackLength_ = vGreat;
        Info<< "    automatic track length specified through"
            << " number of sub cycles : " << nSubCycle_ << nl << endl;
    }
    else
    {
        nSubCycle_ = 1;
        dict.lookup("trackLength") >> trackLength_;
        Info<< "    fixed track length specified : "
            << trackLength_ << nl << endl;
    }

    interpolationScheme_ =
        dict.lookupOrDefault
        (
            "interpolationScheme",
            interpolationCellPoint<scalar>::typeName
        );

    cloudName_ = dict.lookupOrDefault<word>("cloudName", "streamlines");

    meshSearchPtr_.reset(new meshSearch(mesh_));

    sampledSetPtr_ = sampledSet::New
    (
        "seedSampleSet",
        mesh_,
        meshSearchPtr_(),
        dict.subDict("seedSampleSet")
    );

    formatterPtr_ = setWriter::New(dict.lookup("setFormat"), dict);

    return true;
}


Foam::wordList Foam::functionObjects::streamlines::fields() const
{
    wordList allFields(fields_);
    allFields.append(UName_);

    return allFields;
}


bool Foam::functionObjects::streamlines::execute()
{
    return true;
}


bool Foam::functionObjects::streamlines::write()
{
    Info<< type() << " " << name() << " write:" << nl;

    // Create list of available fields
    wordList fieldNames;
    forAll(fields_, fieldi)
    {
        if
        (
            false
            #define FoundTypeField(Type, nullArg) \
              || foundObject<VolField<Type>>(fields_[fieldi])
            FOR_ALL_FIELD_TYPES(FoundTypeField)
            #undef FoundTypeField
        )
        {
            fieldNames.append(fields_[fieldi]);
        }
        else
        {
            cannotFindObject(fields_[fieldi]);
        }
    }

    // Lookup fields and construct interpolators
    #define DeclareTypeInterpolator(Type, nullArg) \
        PtrList<interpolation<Type>> Type##Interp(fieldNames.size());
    FOR_ALL_FIELD_TYPES(DeclareTypeInterpolator);
    #undef DeclareTypeInterpolator
    forAll(fieldNames, fieldi)
    {
        #define ConstructTypeInterpolator(Type, nullArg)                       \
            if (mesh_.foundObject<VolField<Type>>(fieldNames[fieldi]))         \
            {                                                                  \
                Type##Interp.set                                               \
                (                                                              \
                    fieldi,                                                    \
                    interpolation<Type>::New                                   \
                    (                                                          \
                        interpolationScheme_,                                  \
                        mesh_.lookupObject<VolField<Type>>(fieldNames[fieldi]) \
                    )                                                          \
                );                                                             \
            }
        FOR_ALL_FIELD_TYPES(ConstructTypeInterpolator);
        #undef ConstructTypeInterpolator
    }

    // Create a velocity interpolator if it is not already available
    const label UIndex = findIndex(fieldNames, UName_);
    tmpNrc<interpolation<vector>> UInterp(nullptr);
    if (UIndex == -1)
    {
        UInterp =
            tmpNrc<interpolation<vector>>
            (
                interpolation<vector>::New
                (
                    interpolationScheme_,
                    mesh_.lookupObject<volVectorField>(UName_)
                ).ptr()
            );
    }

    // Do tracking to create sampled data
    DynamicField<point> allPositions;
    DynamicField<label> allTracks;
    DynamicField<label> allTrackParts;
    DynamicField<scalar> allAges;
    #define DeclareAllTypes(Type, nullArg) \
        List<DynamicField<Type>> all##Type##s(fieldNames.size());
    FOR_ALL_FIELD_TYPES(DeclareAllTypes);
    #undef DeclareAllTypes
    {
        // Create a cloud and initialise with points from the sampled set
        globalIndex gi(sampledSetPtr_().size());
        streamlinesCloud particles
        (
            mesh_,
            cloudName_,
            IDLList<streamlinesParticle>()
        );
        forAll(sampledSetPtr_(), i)
        {
            particles.addParticle
            (
                new streamlinesParticle
                (
                    mesh_,
                    sampledSetPtr_().positions()[i],
                    sampledSetPtr_().cells()[i],
                    lifeTime_,
                    gi.toGlobal(i)
                )
            );
        }

        // Report the number of successful seeds
        const label nSeeds = returnReduce(particles.size(), sumOp<label>());
        Info << "    Seeded " << nSeeds << " particles" << endl;

        // Create tracking data
        streamlinesParticle::trackingData td
        (
            particles,
            #define TypeInterpolatorParameter(Type, nullArg) \
                Type##Interp,
            FOR_ALL_FIELD_TYPES(TypeInterpolatorParameter)
            #undef TypeInterpolatorParameter
            UIndex != -1 ? vectorInterp[UIndex] : UInterp(),
            trackDirection_ == trackDirection::forward,
            trackOutside_,
            nSubCycle_,
            trackLength_,
            allPositions,
            allTracks,
            allTrackParts,
            allAges
            #define AllTypesParameter(Type, nullArg) \
                , all##Type##s
            FOR_ALL_FIELD_TYPES(AllTypesParameter)
            #undef AllTypesParameter
        );

        // Track
        IDLList<streamlinesParticle> initialParticles;
        if (trackDirection_ == trackDirection::both)
        {
            initialParticles = particles;
        }

        particles.move(particles, td);

        if (trackDirection_ == trackDirection::both)
        {
            particles.IDLList<streamlinesParticle>::operator=(initialParticles);
            td.trackForward_ = !td.trackForward_;
            particles.move(particles, td);
        }
    }

    // Gather data on the master
    if (Pstream::parRun())
    {
        gatherAndFlatten(allPositions);
        gatherAndFlatten(allTracks);
        gatherAndFlatten(allTrackParts);
        gatherAndFlatten(allAges);
        forAll(fieldNames, fieldi)
        {
            #define GatherAndFlattenAllTypes(Type, nullArg) \
                if (Type##Interp.set(fieldi))               \
                {                                           \
                    gatherAndFlatten(all##Type##s[fieldi]); \
                }
            FOR_ALL_FIELD_TYPES(GatherAndFlattenAllTypes);
            #undef GatherAndFlattenAllTypes
        }
    }

    // Report the total number of samples
    Info<< "    Sampled " << allPositions.size() << " locations" << endl;

    // Bin-sort by track and trackPart to build an ordering
    labelList order(allPositions.size());
    if (Pstream::master() && allPositions.size())
    {
        const label nTracks = max(allTracks) + 1;
        const label trackParti0 = min(allTrackParts);
        const label trackParti1 = max(allTrackParts) + 1;

        labelListList trackPartCounts
        (
            nTracks,
            labelList(trackParti1 - trackParti0, 0)
        );
        forAll(allPositions, samplei)
        {
            const label tracki = allTracks[samplei];
            const label trackParti = -trackParti0 + allTrackParts[samplei];
            trackPartCounts[tracki][trackParti] ++;
        }

        label offset = 0;
        labelListList trackPartOffsets
        (
            nTracks,
            labelList(trackParti1 - trackParti0, 0)
        );
        forAll(trackPartOffsets, tracki)
        {
            forAll(trackPartOffsets[tracki], trackParti)
            {
                trackPartOffsets[tracki][trackParti] += offset;
                offset += trackPartCounts[tracki][trackParti];
            }
        }

        forAll(trackPartCounts, tracki)
        {
            trackPartCounts[tracki] = 0;
        }

        forAll(allPositions, samplei)
        {
            const label tracki = allTracks[samplei];
            const label trackParti = -trackParti0 + allTrackParts[samplei];

            order[samplei] =
                trackPartOffsets[tracki][trackParti]
              + trackPartCounts[tracki][trackParti];

            trackPartCounts[tracki][trackParti] ++;
        }
    }

    //auto reportTrackParts = [&]()
    //{
    //    Info<< nl;
    //    forAll(allPositions, samplei)
    //    {
    //        if
    //        (
    //            samplei == 0
    //         || allTracks[samplei] != allTracks[samplei - 1]
    //         || allTrackParts[samplei] != allTrackParts[samplei - 1]
    //        )
    //        {
    //            Info<< "track #" << allTracks[samplei]
    //                << " part #" << allTrackParts[samplei]
    //                << " from i=" << samplei << " to ";
    //        }
    //        if
    //        (
    //            samplei == allPositions.size() - 1
    //         || allTracks[samplei + 1] != allTracks[samplei]
    //         || allTrackParts[samplei + 1] != allTrackParts[samplei]
    //        )
    //        {
    //            Info<< "i=" << samplei << nl;
    //        }
    //    }
    //};

    //reportTrackParts();

    // Reorder
    if (Pstream::master())
    {
        allPositions.rmap(allPositions, order);
        allTracks.rmap(allTracks, order);
        allTrackParts.rmap(allTrackParts, order);
        allAges.rmap(allAges, order);
        forAll(fieldNames, fieldi)
        {
            #define RMapAllTypes(Type, nullArg)                             \
                if (Type##Interp.set(fieldi))                               \
                {                                                           \
                    all##Type##s[fieldi].rmap(all##Type##s[fieldi], order); \
                }
            FOR_ALL_FIELD_TYPES(RMapAllTypes);
            #undef RMapAllTypes
        }
    }

    //reportTrackParts();

    // Relabel tracks and track parts into track labels only, and join the
    // forward and backward track parts that are connected to the seed
    if (Pstream::master())
    {
        label samplei = 0, tracki = 0;
        forAll(allPositions, samplej)
        {
            const label trackj = allTracks[samplej];
            const label trackPartj = allTrackParts[samplej];

            allPositions[samplei] = allPositions[samplej];
            allTracks[samplei] = tracki;
            allTrackParts[samplei] = 0;
            allAges[samplei] = allAges[samplej];
            forAll(fieldNames, fieldi)
            {
                #define ShuffleUpAllTypes(Type, nullArg)   \
                    if (Type##Interp.set(fieldi))          \
                    {                                      \
                        all##Type##s[fieldi][samplei] =    \
                            all##Type##s[fieldi][samplej]; \
                    }
                FOR_ALL_FIELD_TYPES(ShuffleUpAllTypes);
                #undef ShuffleUpAllTypes
            }

            const bool joinNewParts =
                samplej != allPositions.size() - 1
             && trackPartj == -1
             && allTrackParts[samplej + 1] == 0;

            if (!joinNewParts) samplei ++;

            const bool newPart =
                samplej == allPositions.size() - 1
             || trackj != allTracks[samplej + 1]
             || trackPartj != allTrackParts[samplej + 1];

            if (!joinNewParts && newPart) tracki ++;
        }

        allPositions.resize(samplei);
        allTracks.resize(samplei);
        allTrackParts.resize(samplei);
        allAges.resize(samplei);
        forAll(fieldNames, fieldi)
        {
            #define ResizeAllTypes(Type, nullArg)         \
                if (Type##Interp.set(fieldi))             \
                {                                         \
                    all##Type##s[fieldi].resize(samplei); \
                }
            FOR_ALL_FIELD_TYPES(ResizeAllTypes);
            #undef ResizeAllTypes
        }
    }

    //reportTrackParts();

    // Write
    if (Pstream::master() && allPositions.size())
    {
        // Make output directory
        const fileName outputPath =
            time_.globalPath()
           /writeFile::outputPrefix
           /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
           /name()
           /time_.name();
        mkDir(outputPath);

        // Pass data to the formatter to write
        const label nValueSets = fieldNames.size() + writeAge_;
        wordList valueSetNames(nValueSets);
        #define DeclareTypeValueSets(Type, nullArg) \
            UPtrList<const Field<Type>> Type##ValueSets(nValueSets);
        FOR_ALL_FIELD_TYPES(DeclareTypeValueSets);
        #undef DeclareTypeValueSets
        if (writeAge_)
        {
            valueSetNames[0] = "age";
            scalarValueSets.set(0, &allAges);
        }
        forAll(fieldNames, fieldi)
        {
            valueSetNames[fieldi + writeAge_] = fieldNames[fieldi];

            #define SetTypeValueSetPtr(Type, nullArg) \
                if (Type##Interp.set(fieldi))         \
                {                                     \
                    Type##ValueSets.set               \
                    (                                 \
                        fieldi + writeAge_,           \
                        &all##Type##s[fieldi]         \
                    );                                \
                }
            FOR_ALL_FIELD_TYPES(SetTypeValueSetPtr);
            #undef SetTypeValueSetPtr
        }
        formatterPtr_->write
        (
            outputPath,
            "tracks",
            coordSet(allTracks, word::null, allPositions),
            valueSetNames
            #define TypeValueSetsParameter(Type, nullArg) , Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
            #undef TypeValueSetsParameter
        );
    }

    return true;
}


void Foam::functionObjects::streamlines::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        // Moving mesh affects the search tree
        read(dict_);
    }
}


void Foam::functionObjects::streamlines::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        read(dict_);
    }
}


void Foam::functionObjects::streamlines::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        read(dict_);
    }
}


void Foam::functionObjects::streamlines::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        read(dict_);
    }
}


// ************************************************************************* //
