/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "vtkPV3Foam.H"

// OpenFOAM includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "Cloud.H"
#include "vtkPV3FoamReader.h"

// local headers
#include "vtkPV3FoamAddToSelection.H"
#include "vtkPV3FoamUpdateInfoFields.H"

// VTK includes
#include "vtkDataArraySelection.h"


// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
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
            PtrList<entry>(readStream("regIOobject"))
        {
            close();
        }

   // Member functions

        bool writeData(Ostream&) const
        {
            notImplemented("zonesEntries::writeData(Ostream&) const");
            return false;
        }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType>
Foam::wordList Foam::vtkPV3Foam::getZoneNames
(
    const ZoneMesh<ZoneType, polyMesh>& zmesh
) const
{
    wordList names(zmesh.size());
    label nZone = 0;

    forAll(zmesh, zoneI)
    {
        if (zmesh[zoneI].size())
        {
            names[nZone++] = zmesh[zoneI].name();
        }
    }
    names.setSize(nZone);

    return names;
}


Foam::wordList Foam::vtkPV3Foam::getZoneNames(const word& zoneType) const
{
    wordList names;

    // mesh not loaded - read from file
    IOobject ioObj
    (
        zoneType,
        dbPtr_().findInstance
        (
            meshDir_,
            zoneType,
            IOobject::READ_IF_PRESENT
        ),
        meshDir_,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.headerOk())
    {
        zonesEntries zones(ioObj);

        names.setSize(zones.size());
        forAll(zones, zoneI)
        {
            names[zoneI] = zones[zoneI].keyword();
        }
    }

    return names;
}


void Foam::vtkPV3Foam::updateInfoInternalMesh
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoInternalMesh" << endl;
    }

    // Determine mesh parts (internalMesh, patches...)
    //- Add internal mesh as first entry
    arrayRangeVolume_.reset(arraySelection->GetNumberOfArrays());
    arraySelection->AddArray
    (
        "internalMesh"
    );
    arrayRangeVolume_ += 1;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoInternalMesh" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoLagrangian
(
    vtkDataArraySelection* arraySelection
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoLagrangian" << nl
            << "    " << dbPtr_->timePath()/cloud::prefix << endl;
    }


    // use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    // Search for list of lagrangian objects for this time
    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/lagrangianPrefix, fileName::DIRECTORY)
    );

    arrayRangeLagrangian_.reset(arraySelection->GetNumberOfArrays());

    int nClouds = 0;
    forAll(cloudDirs, cloudI)
    {
        // Add cloud to GUI list
        arraySelection->AddArray
        (
            (cloudDirs[cloudI] + " - lagrangian").c_str()
        );

        ++nClouds;
    }
    arrayRangeLagrangian_ += nClouds;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoLagrangian" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoPatches
(
    vtkDataArraySelection* arraySelection,
    stringList& enabledEntries
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoPatches"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }


    HashSet<string> enabledEntriesSet(enabledEntries);

    arrayRangePatches_.reset(arraySelection->GetNumberOfArrays());

    int nPatches = 0;
    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        const HashTable<labelList, word>& groups = patches.groupPatchIDs();
        const wordList allPatchNames = patches.names();

        // Add patch groups
        // ~~~~~~~~~~~~~~~~

        for
        (
            HashTable<labelList, word>::const_iterator iter = groups.begin();
            iter != groups.end();
            ++iter
        )
        {
            const word& groupName = iter.key();
            const labelList& patchIDs = iter();

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
                                string vtkPatchName = pp.name() + " - patch";
                                enabledEntriesSet.insert(vtkPatchName);
                            }
                        }
                    }
                }
            }
        }


        // Add patches
        // ~~~~~~~~~~~

        if (!reader_->GetShowGroupsOnly())
        {
            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (pp.size())
                {
                    // Add patch to GUI list
                    arraySelection->AddArray
                    (
                        (pp.name() + " - patch").c_str()
                    );

                    ++nPatches;
                }
            }
        }
    }
    else
    {
        // mesh not loaded - read from file
        // but this could fail if we've supplied a bad region name
        IOobject ioObj
        (
            "boundary",
            dbPtr_().findInstance
            (
                meshDir_,
                "boundary",
                IOobject::READ_IF_PRESENT
            ),
            meshDir_,
            dbPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        // this should only ever fail if the mesh region doesn't exist
        if (ioObj.headerOk())
        {
            polyBoundaryMeshEntries patchEntries(ioObj);


            // Read patches and determine sizes
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            wordList names(patchEntries.size());
            labelList sizes(patchEntries.size());

            forAll(patchEntries, patchI)
            {
                const dictionary& patchDict = patchEntries[patchI].dict();

                sizes[patchI] = readLabel(patchDict.lookup("nFaces"));
                names[patchI] = patchEntries[patchI].keyword();
            }


            // Add (non-zero) patch groups to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            HashTable<labelList, word> groups(patchEntries.size());

            forAll(patchEntries, patchI)
            {
                const dictionary& patchDict = patchEntries[patchI].dict();

                wordList groupNames;
                patchDict.readIfPresent("inGroups", groupNames);

                forAll(groupNames, groupI)
                {
                    HashTable<labelList, word>::iterator iter = groups.find
                    (
                        groupNames[groupI]
                    );
                    if (iter != groups.end())
                    {
                        iter().append(patchI);
                    }
                    else
                    {
                        groups.insert(groupNames[groupI], labelList(1, patchI));
                    }
                }
            }

            for
            (
                HashTable<labelList, word>::const_iterator iter =
                    groups.begin();
                iter != groups.end();
                ++iter
            )
            {
                const word& groupName = iter.key();
                const labelList& patchIDs = iter();

                label nFaces = 0;
                forAll(patchIDs, i)
                {
                    nFaces += sizes[patchIDs[i]];
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
                                    string vtkPatchName =
                                        names[patchIDs[i]] + " - patch";
                                    enabledEntriesSet.insert(vtkPatchName);
                                }
                            }
                        }
                    }
                }
            }


            // Add (non-zero) patches to the list of mesh parts
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (!reader_->GetShowGroupsOnly())
            {
                forAll(names, patchI)
                {
                    // Valid patch if nFace > 0 - add patch to GUI list
                    if (sizes[patchI])
                    {
                        arraySelection->AddArray
                        (
                            (names[patchI] + " - patch").c_str()
                        );

                        ++nPatches;
                    }
                }
            }
        }
    }
    arrayRangePatches_ += nPatches;

    // Update enabled entries in case of group selection
    enabledEntries = enabledEntriesSet.toc();

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoPatches" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoZones
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeZones())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoZones"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    wordList namesLst;

    //
    // cellZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->cellZones());
    }
    else
    {
        namesLst = getZoneNames("cellZones");
    }

    arrayRangeCellZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            (namesLst[elemI] + " - cellZone").c_str()
        );
    }
    arrayRangeCellZones_ += namesLst.size();


    //
    // faceZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->faceZones());
    }
    else
    {
        namesLst = getZoneNames("faceZones");
    }

    arrayRangeFaceZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            (namesLst[elemI] + " - faceZone").c_str()
        );
    }
    arrayRangeFaceZones_ += namesLst.size();


    //
    // pointZones information
    // ~~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = getZoneNames(meshPtr_->pointZones());
    }
    else
    {
        namesLst = getZoneNames("pointZones");
    }

    arrayRangePointZones_.reset(arraySelection->GetNumberOfArrays());
    forAll(namesLst, elemI)
    {
        arraySelection->AddArray
        (
            (namesLst[elemI] + " - pointZone").c_str()
        );
    }
    arrayRangePointZones_ += namesLst.size();

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoZones" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoSets
(
    vtkDataArraySelection* arraySelection
)
{
    if (!reader_->GetIncludeSets())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoSets" << endl;
    }

    // Add names of sets. Search for last time directory with a sets
    // subdirectory. Take care not to search beyond the last mesh.

    word facesInstance = dbPtr_().findInstance
    (
        meshDir_,
        "faces",
        IOobject::READ_IF_PRESENT
    );

    word setsInstance = dbPtr_().findInstance
    (
        meshDir_/"sets",
        word::null,
        IOobject::READ_IF_PRESENT,
        facesInstance
    );

    IOobjectList objects(dbPtr_(), setsInstance, meshDir_/"sets");

    if (debug)
    {
        Info<< "     Foam::vtkPV3Foam::updateInfoSets read "
            << objects.names() << " from " << setsInstance << endl;
    }


    arrayRangeCellSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangeCellSets_ += addToSelection<cellSet>
    (
        arraySelection,
        objects,
        " - cellSet"
    );

    arrayRangeFaceSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangeFaceSets_ += addToSelection<faceSet>
    (
        arraySelection,
        objects,
        " - faceSet"
    );

    arrayRangePointSets_.reset(arraySelection->GetNumberOfArrays());
    arrayRangePointSets_ += addToSelection<pointSet>
    (
        arraySelection,
        objects,
        " - pointSet"
    );

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(arraySelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoSets" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoLagrangianFields()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoLagrangianFields"
            << endl;
    }

    vtkDataArraySelection* fieldSelection =
        reader_->GetLagrangianFieldSelection();

    // preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(fieldSelection);
    fieldSelection->RemoveAllArrays();

    // TODO - currently only get fields from ONE cloud
    // have to decide if the second set of fields get mixed in
    // or dealt with separately

    const arrayRange& range = arrayRangeLagrangian_;
    if (range.empty())
    {
        return;
    }

    int partId = range.start();
    word cloudName = getPartName(partId);

    // use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().timeName(),
        lagrangianPrefix/cloudName
    );

    addToSelection<IOField<label> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<scalar> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<vector> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<sphericalTensor> >
    (
        fieldSelection,

        objects
    );
    addToSelection<IOField<symmTensor> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<tensor> >
    (
        fieldSelection,
        objects
    );

    // restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateInfoLagrangianFields - "
            << "lagrangian objects.size() = " << objects.size() << endl;
    }
}


// ************************************************************************* //
