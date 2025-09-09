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

\*---------------------------------------------------------------------------*/

#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"
#include "vtkPVFoamAddToSelection.H"

// OpenFOAM includes
#include "processorRunTimes.H"
#include "domainDecomposition.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "cloud.H"
#include "pointMesh.H"
#include "LagrangianMesh.H"

// VTK includes
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
template<class ZonesType>
class zonesEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    // Constructors

        explicit zonesEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream(ZonesType::typeName))
        {
            close();
        }

   // Member functions

        bool writeData(Ostream&) const
        {
            NotImplemented;
            return false;
        }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZonesType>
Foam::wordList Foam::vtkPVFoam::getZoneNames(const ZonesType& zones) const
{
    wordList names(zones.size());
    label nZone = 0;

    forAll(zones, zoneI)
    {
        if (zones[zoneI].size())
        {
            names[nZone++] = zones[zoneI].name();
        }
    }
    names.setSize(nZone);

    return names;
}


template<class ZonesType>
Foam::wordList Foam::vtkPVFoam::getZoneNames(const word& zonesName) const
{
    const Time& runTime =
        reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time()
      : procDbsPtr_->completeTime();

    // Mesh not loaded - read from file
    typeIOobject<ZonesType> ioObj
    (
        zonesName,
        runTime.findInstance
        (
            meshDir_,
            zonesName,
            IOobject::READ_IF_PRESENT
        ),
        meshDir_,
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.headerOk())
    {
        const zonesEntries<ZonesType> zones(ioObj);

        wordList names(zones.size());

        forAll(zones, zoneI)
        {
            names[zoneI] = zones[zoneI].keyword();
        }

        return names;
    }
    else
    {
        return wordList();
    }
}


void Foam::vtkPVFoam::updateInfoInternalMesh
(
    vtkDataArraySelection* arraySelection
)
{
    // Determine mesh parts (internalMesh, patches...)
    // - Add internal mesh as first entry
    arrayRangeVolume_.reset(arraySelection->GetNumberOfArrays());
    arraySelection->AddArray
    (
        "internalMesh"
    );
    arrayRangeVolume_ += 1;
}


void Foam::vtkPVFoam::updateInfolagrangian
(
    vtkDataArraySelection* arraySelection
)
{
    UPtrList<const Time> runTimes;
    if (reader_->GetDecomposedCase())
    {
        runTimes.setSize(procDbsPtr_->nProcs());
        forAll(procDbsPtr_->procTimes(), proci)
        {
            runTimes.set(proci, &procDbsPtr_->procTimes()[proci]);
        }
    }
    else
    {
        runTimes.setSize(1);
        runTimes.set(0, &procDbsPtr_->completeTime());
    }

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    const fileName lagrangianPrefix =
        meshRegion_ == polyMesh::defaultRegion
      ? fileName(lagrangian::cloud::prefix)
      : meshRegion_/lagrangian::cloud::prefix;

    arrayRangelagrangian_.reset(arraySelection->GetNumberOfArrays());

    // Generate a list of lagrangian clouds across all times
    HashSet<fileName> lagrangianDirs;

    // Get times list. Flush first to force refresh.
    fileHandler().flush();
    forAll(runTimes, runTimei)
    {
        const instantList times = runTimes[runTimei].times();

        forAll(times, timei)
        {
            lagrangianDirs +=
                fileHandler().readDir
                (
                    fileHandler().filePath
                    (
                        runTimes[runTimei].path()
                       /times[timei].name()
                       /lagrangianPrefix
                    ),
                    fileType::directory
                );
        }
    }

    forAllConstIter(HashSet<fileName>, lagrangianDirs, cloudIter)
    {
        arraySelection->AddArray
        (
            (cloudIter.key() + " - lagrangian").c_str()
        );
    }

    arrayRangelagrangian_ += lagrangianDirs.size();
}


void Foam::vtkPVFoam::updateInfoLagrangian
(
    vtkDataArraySelection* arraySelection
)
{
    const Time& runTime =
        reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time()
      : procDbsPtr_->completeTime();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    const fileName LagrangianPrefix =
        meshRegion_ == polyMesh::defaultRegion
      ? fileName(LagrangianMesh::prefix)
      : meshRegion_/LagrangianMesh::prefix;

    arrayRangeLagrangian_.reset(arraySelection->GetNumberOfArrays());

    // Generate a list of Lagrangian meshes across all times
    HashSet<fileName> LagrangianDirs;

    // Get times list. Flush first to force refresh.
    fileHandler().flush();
    instantList times = runTime.times();
    forAll(times, timei)
    {
        LagrangianDirs +=
            fileHandler().readDir
            (
                fileHandler().filePath
                (
                    runTime.path()
                   /times[timei].name()
                   /LagrangianPrefix
                ),
                fileType::directory
            );
    }

    forAllConstIter(HashSet<fileName>, LagrangianDirs, LagrangianIter)
    {
        arraySelection->AddArray
        (
            (LagrangianIter.key() + " - Lagrangian").c_str()
        );
    }

    arrayRangeLagrangian_ += LagrangianDirs.size();
}


void Foam::vtkPVFoam::updateInfoPatches
(
    vtkDataArraySelection* arraySelection,
    stringList& enabledEntries,
    const bool first
)
{
    HashSet<string> enabledEntriesSet(enabledEntries);

    arrayRangePatches_.reset(arraySelection->GetNumberOfArrays());

    int nPatches = 0;
    if (procMeshesPtr_.valid())
    {
        const fvMesh& mesh = procMeshesPtr_->completeMesh();

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const HashTable<labelList, word>& groups = patches.groupPatchIndices();
        const wordList allPatchNames = patches.names();

        // Add patch groups
        for
        (
            HashTable<labelList, word>::const_iterator iter = groups.begin();
            iter != groups.end();
            ++iter
        )
        {
            const word& groupName = iter.key();
            const labelList& patchIDs = iter();

            // Count the number of faces in this group
            label nFaces = 0;
            forAll(patchIDs, i)
            {
                nFaces += patches[patchIDs[i]].size();
            }

            // Valid patch if nFace > 0 - add patch to GUI list
            if (nFaces)
            {
                string vtkGrpName = groupName + " - group";
                arraySelection->AddArray(vtkGrpName.c_str());

                ++nPatches;

                if (enabledEntriesSet.found(vtkGrpName))
                {
                    if (!reader_->GetShowGroupsOnly())
                    {
                        enabledEntriesSet.erase(vtkGrpName);
                        forAll(patchIDs, i)
                        {
                            const polyPatch& pp = patches[patchIDs[i]];

                            if (pp.size())
                            {
                                enabledEntriesSet.insert
                                (
                                    pp.name() + " - " + pp.type()
                                );
                            }
                        }
                    }
                }
            }
        }

        // Add patches
        if (!reader_->GetShowGroupsOnly())
        {
            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi];

                if (pp.size())
                {
                    arraySelection->AddArray
                    (
                        (pp.name() + " - " + pp.type()).c_str()
                    );

                    ++nPatches;
                }
            }
        }
    }
    else
    {
        const bool decomposed = reader_->GetDecomposedCase();

        const Time& runTime =
            decomposed
          ? procDbsPtr_->proc0Time()
          : procDbsPtr_->completeTime();

        // Mesh not loaded - read from file
        // but this could fail if we've supplied a bad region name
        typeIOobject<polyBoundaryMesh> ioObj
        (
            "boundary",
            runTime.findInstance
            (
                meshDir_,
                "boundary",
                IOobject::READ_IF_PRESENT
            ),
            meshDir_,
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        // This should only ever fail if the mesh region doesn't exist
        if (ioObj.headerOk())
        {
            polyBoundaryMeshEntries patchEntries(ioObj);

            // Read patches and determine sizes
            wordList names(patchEntries.size());
            labelList sizes(patchEntries.size());
            boolList isProc(patchEntries.size());
            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();

                names[patchi] = patchEntries[patchi].keyword();
                sizes[patchi] = patchDict.lookup<label>("nFaces");
                isProc[patchi] = patchDict.found("myProcNo");
            }

            // Build a map from the group name to a list of the indices of the
            // patches in the group
            HashTable<labelList> groups(patchEntries.size());
            forAll(patchEntries, patchi)
            {
                const dictionary& patchDict = patchEntries[patchi].dict();

                const wordList groupNames =
                    patchDict.lookupOrDefault("inGroups", wordList());

                forAll(groupNames, groupI)
                {
                    HashTable<labelList, word>::iterator iter =
                        groups.find(groupNames[groupI]);

                    if (iter != groups.end())
                    {
                        iter().append(patchi);
                    }
                    else
                    {
                        groups.insert(groupNames[groupI], labelList(1, patchi));
                    }
                }
            }

            // Add (non-zero) patch groups to the list of mesh parts
            forAllConstIter(HashTable<labelList>, groups, iter)
            {
                const word& groupName = iter.key();
                const labelList& patchIDs = iter();

                label nFaces = 0;
                forAll(patchIDs, i)
                {
                    if (!isProc[patchIDs[i]])
                    {
                        nFaces += sizes[patchIDs[i]];
                    }
                }

                // Valid patch if nFace > 0 - add patch to GUI list
                if (nFaces)
                {
                    string vtkGrpName = groupName + " - group";
                    arraySelection->AddArray(vtkGrpName.c_str());

                    ++nPatches;

                    if (enabledEntriesSet.found(vtkGrpName))
                    {
                        if (!reader_->GetShowGroupsOnly())
                        {
                            enabledEntriesSet.erase(vtkGrpName);

                            forAll(patchIDs, i)
                            {
                                if (sizes[patchIDs[i]])
                                {
                                    enabledEntriesSet.insert
                                    (
                                        names[patchIDs[i]]
                                      + " - "
                                      + patchEntries[patchIDs[i]]
                                       .dict()
                                       .lookup<word>("type")
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // Add (non-zero) patches to the list of mesh parts
            if (!reader_->GetShowGroupsOnly())
            {
                const wordReList defaultPatchTypes
                (
                    configDict_.lookupOrDefault
                    (
                        "defaultPatchTypes",
                        wordReList(wordList({"patch", "wall"}))
                    )
                );

                forAll(names, patchi)
                {
                    // Valid patch if nFace > 0 - add patch to GUI list
                    if (sizes[patchi])
                    {
                        const word patchType
                        (
                            patchEntries[patchi].dict().lookup("type")
                        );

                        const string vtkPatchName
                        (
                            names[patchi] + " - " + patchType
                        );

                        arraySelection->AddArray(vtkPatchName.c_str());

                        if (first && findStrings(defaultPatchTypes, patchType))
                        {
                            enabledEntriesSet.insert(vtkPatchName);
                        }

                        ++nPatches;
                    }
                }
            }
        }
    }

    arrayRangePatches_ += nPatches;

    // Update enabled entries in case of group selection
    enabledEntries = enabledEntriesSet.toc();
}


void Foam::vtkPVFoam::updateInfoZones
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeZones()) return;

    const fvMesh& mesh =
       !procMeshesPtr_.valid()
      ? NullObjectRef<fvMesh>()
      : reader_->GetDecomposedCase()
      ? procMeshesPtr_->haveProcs()
      ? procMeshesPtr_->procMeshes().first()
      : NullObjectRef<fvMesh>()
      : procMeshesPtr_->completeMesh();

    // cellZones information
    const wordList cellZoneNames =
        notNull(mesh)
      ? getZoneNames(mesh.cellZones())
      : getZoneNames<cellZoneList>("cellZones");

    arrayRangeCellZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(cellZoneNames, i)
    {
        arraySelection->AddArray
        (
            (cellZoneNames[i] + " - cellZone").c_str()
        );
    }
    arrayRangeCellZones_ += cellZoneNames.size();

    // faceZones information
    const wordList faceZoneNames =
        notNull(mesh)
      ? getZoneNames(mesh.faceZones())
      : getZoneNames<faceZoneList>("faceZones");

    arrayRangeFaceZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(faceZoneNames, i)
    {
        arraySelection->AddArray
        (
            (faceZoneNames[i] + " - faceZone").c_str()
        );
    }
    arrayRangeFaceZones_ += faceZoneNames.size();

    // pointZones information
    const wordList pointZoneNames =
        notNull(mesh)
      ? getZoneNames(mesh.pointZones())
      : getZoneNames<pointZoneList>("pointZones");

    arrayRangePointZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(pointZoneNames, i)
    {
        arraySelection->AddArray
        (
            (pointZoneNames[i] + " - pointZone").c_str()
        );
    }
    arrayRangePointZones_ += pointZoneNames.size();
}


void Foam::vtkPVFoam::updateInfoSets
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeSets()) return;

    const Time& runTime =
        reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time()
      : procDbsPtr_->completeTime();

    // Add names of sets. Search for last time directory with a sets
    // subdirectory. Take care not to search beyond the last mesh.

    const word facesInstance = runTime.findInstance
    (
        meshDir_,
        "faces",
        IOobject::READ_IF_PRESENT
    );

    const word setsInstance = runTime.findInstance
    (
        meshDir_/"sets",
        word::null,
        IOobject::READ_IF_PRESENT,
        facesInstance
    );

    const IOobjectList objects(runTime, setsInstance, meshDir_/"sets");

    arrayRangeCellSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangeCellSets_ +=
        addToSelection<cellSet>(arraySelection, objects, " - cellSet");

    arrayRangeFaceSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangeFaceSets_ +=
        addToSelection<faceSet>(arraySelection, objects, " - faceSet");

    arrayRangePointSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangePointSets_ +=
        addToSelection<pointSet>(arraySelection, objects, " - pointSet");
}


void Foam::vtkPVFoam::updateInfoFields()
{
    DebugInFunction;

    vtkDataArraySelection* fieldSelection = reader_->GetFieldSelection();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    word regionPrefix;
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        regionPrefix = meshRegion_;
    }

    const Time& runTime =
        reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time()
      : procDbsPtr_->completeTime();

    const instantList times = runTime.times();

    stringList enabledEntries;
    if (fieldSelection->GetNumberOfArrays() == 0 && !procMeshesPtr_.valid())
    {
        const wordReList defaultFieldRes
        (
            configDict_.lookupOrDefault
            (
                "defaultFields",
                wordReList(wordList{"p", "U"})
            )
        );

        wordHashSet objectNameSet;
        forAll(times, timei)
        {
            // Search for list of objects for this time and mesh region
            IOobjectList objects(runTime, times[timei].name(), regionPrefix);

            forAllConstIter(IOobjectList, objects, iter)
            {
                objectNameSet.insert(iter.key());
            }
        }

        const wordList objectNames(objectNameSet.toc());
        const labelList defaultFields
        (
            findStrings(defaultFieldRes, objectNames)
        );

        enabledEntries.setSize(defaultFields.size());
        forAll(defaultFields, i)
        {
            enabledEntries[i] = objectNames[defaultFields[i]];
        }
    }
    else
    {
        // Preserve the enabled selections
        enabledEntries = getSelectedArrayEntries(fieldSelection, false);
    }

    fieldSelection->RemoveAllArrays();

    forAll(times, timei)
    {
        // Search for list of objects for this time and mesh region
        IOobjectList objects(runTime, times[timei].name(), regionPrefix);

        addFieldsToSelection<volMesh>(fieldSelection, objects);
        addInternalFieldsToSelection<volMesh>(fieldSelection, objects);
        addFieldsToSelection<surfaceMesh>(fieldSelection, objects);
        addFieldsToSelection<pointMesh>(fieldSelection, objects);
    }

    // Restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);

    if (debug) getSelectedArrayEntries(fieldSelection);
}


void Foam::vtkPVFoam::updateInfolagrangianFields()
{
    DebugInFunction;

    UPtrList<const Time> runTimes;
    if (reader_->GetDecomposedCase())
    {
        runTimes.setSize(procDbsPtr_->nProcs());
        forAll(procDbsPtr_->procTimes(), proci)
        {
            runTimes.set(proci, &procDbsPtr_->procTimes()[proci]);
        }
    }
    else
    {
        runTimes.setSize(1);
        runTimes.set(0, &procDbsPtr_->completeTime());
    }

    vtkDataArraySelection* fieldSelection =
        reader_->GetlagrangianFieldSelection();

    // Preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(fieldSelection, false);
    fieldSelection->RemoveAllArrays();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(lagrangian::cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/lagrangian::cloud::prefix;
    }

    // Add the available fields from all clouds and all time directories.
    // Differing sets of fields from multiple clouds get combined into a single
    // set. ParaView will display "(partial)" after field names that only apply
    // to some of the clouds.
    const arrayRange& range = arrayRangelagrangian_;
    fileHandler().flush();
    for (label partId = range.start(); partId < range.end(); ++ partId)
    {
        forAll(runTimes, runTimei)
        {
            const instantList times = runTimes[runTimei].times();

            forAll(times, timei)
            {
                IOobjectList objects
                (
                    runTimes[runTimei],
                    times[timei].name(),
                    lagrangianPrefix/getPartName(partId)
                );

                #define ADD_TO_SELECTION(Type, nullArg) \
                    addToSelection<IOField<Type>>(fieldSelection, objects);
                ADD_TO_SELECTION(label, );
                FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION)
                #undef ADD_TO_SELECTION
            }
        }
    }

    // Restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);

    if (debug) getSelectedArrayEntries(fieldSelection);
}


void Foam::vtkPVFoam::updateInfoLagrangianFields()
{
    DebugInFunction;

    const Time& runTime =
        reader_->GetDecomposedCase()
      ? procDbsPtr_->proc0Time()
      : procDbsPtr_->completeTime();

    const instantList times = runTime.times();

    vtkDataArraySelection* fieldSelection =
        reader_->GetLagrangianFieldSelection();

    // Preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(fieldSelection, false);
    fieldSelection->RemoveAllArrays();

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    const fileName LagrangianPrefix =
        meshRegion_ == polyMesh::defaultRegion
      ? fileName(LagrangianMesh::prefix)
      : meshRegion_/LagrangianMesh::prefix;

    // Add the available fields from all clouds and all time directories.
    // Differing sets of fields from multiple clouds get combined into a single
    // set. ParaView will display "(partial)" after field names that only apply
    // to some of the clouds.
    const arrayRange& range = arrayRangeLagrangian_;

    fileHandler().flush();

    for (label partId = range.start(); partId < range.end(); ++ partId)
    {
        forAll(times, timei)
        {
            IOobjectList objects
            (
                runTime,
                times[timei].name(),
                LagrangianPrefix/getPartName(partId)
            );

            #define ADD_TO_SELECTION(Type, GeoField) \
                addToSelection<GeoField<Type>>(fieldSelection, objects);
            ADD_TO_SELECTION(label, LagrangianField);
            FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION, LagrangianField);
            ADD_TO_SELECTION(label, LagrangianInternalField);
            FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION, LagrangianInternalField);
            #undef ADD_TO_SELECTION
        }
    }

    // Restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);

    if (debug) getSelectedArrayEntries(fieldSelection);
}


// ************************************************************************* //
