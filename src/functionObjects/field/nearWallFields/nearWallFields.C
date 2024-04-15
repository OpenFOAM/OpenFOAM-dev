/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "nearWallFields.H"
#include "wordReList.H"
#include "findCellParticle.H"
#include "OBJstream.H"
#include "globalIndex.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(nearWallFields, 0);
    addToRunTimeSelectionTable(functionObject, nearWallFields, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::nearWallFields::calcAddressing()
{
    // Count number of faces
    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        nPatchFaces += mesh_.boundary()[patchi].size();
    }

    // Global indexing
    globalIndex globalWalls(nPatchFaces);

    DebugInFunction << "nPatchFaces: " << globalWalls.size() << endl;

    // Construct cloud
    Cloud<findCellParticle> cloud
    (
        mesh_,
        cloud::defaultName,
        IDLList<findCellParticle>()
    );

    // Add particles to track to sample locations
    nPatchFaces = 0;

    label nLocateBoundaryHits = 0;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& patch = mesh_.boundary()[patchi];

        vectorField nf(patch.nf());

        forAll(patch, patchFacei)
        {
            cloud.addParticle
            (
                new findCellParticle
                (
                    mesh_,
                    patch.Cf()[patchFacei],
                    patch.faceCells()[patchFacei],
                    nLocateBoundaryHits,
                  - distance_*nf[patchFacei],
                    globalWalls.toGlobal(nPatchFaces) // passive data
                )
            );

            nPatchFaces++;
        }
    }



    if (debug)
    {
        // Dump particles
        OBJstream str
        (
            time_.path()
           /"wantedTracks_" + time_.name() + ".obj"
        );
        InfoInFunction << "Dumping tracks to " << str.name() << endl;

        forAllConstIter(Cloud<findCellParticle>, cloud, iter)
        {
            const vector p = iter().position(mesh_);
            str.write(linePointRef(p, p + iter().displacement()));
        }
    }



    // Per cell: empty or global wall index and end location
    cellToWalls_.setSize(mesh_.nCells());
    cellToSamples_.setSize(mesh_.nCells());

    // Database to pass into findCellParticle::move
    findCellParticle::trackingData td(cloud, cellToWalls_, cellToSamples_);

    // Debug: collect start points
    pointField start;
    if (debug)
    {
        start.setSize(nPatchFaces);
        nPatchFaces = 0;
        forAllConstIter(Cloud<findCellParticle>, cloud, iter)
        {
            const findCellParticle& tp = iter();
            start[nPatchFaces++] = tp.position(mesh_);
        }
    }


    // Track
    cloud.move(cloud, td);


    // Rework cell-to-globalpatchface into a map
    List<Map<label>> compactMap;
    getPatchDataMapPtr_.reset
    (
        new distributionMap
        (
            globalWalls,
            cellToWalls_,
            compactMap
        )
    );


    // Debug: dump resulting tracks
    if (debug)
    {
        getPatchDataMapPtr_().distribute(start);
        {
            OBJstream str
            (
                time_.path()
               /"obtainedTracks_" + time_.name() + ".obj"
            );
            InfoInFunction << "Dumping obtained to " << str.name() << endl;

            forAll(cellToWalls_, celli)
            {
                const List<point>& ends = cellToSamples_[celli];
                const labelList& cData = cellToWalls_[celli];
                forAll(cData, i)
                {
                    str.write(linePointRef(ends[i], start[cData[i]]));
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::nearWallFields::nearWallFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::nearWallFields::~nearWallFields()
{
    DebugInFunction << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::nearWallFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("fields") >> fieldSet_;
    patchSet_ = patchSet(dict);
    distance_ = dict.lookup<scalar>("distance");


    // Clear out any previously loaded fields
    vsf_.clear();
    vvf_.clear();
    vSpheretf_.clear();
    vSymmtf_.clear();
    vtf_.clear();
    fieldMap_.clear();
    reverseFieldMap_.clear();


    // Generate fields with mappedField boundary condition

    // Convert field to map
    fieldMap_.resize(2*fieldSet_.size());
    reverseFieldMap_.resize(2*fieldSet_.size());
    forAll(fieldSet_, seti)
    {
        const word& fieldName = fieldSet_[seti].first();
        const word& sampleFieldName = fieldSet_[seti].second();

        fieldMap_.insert(fieldName, sampleFieldName);
        reverseFieldMap_.insert(sampleFieldName, fieldName);
    }

    Log  << type() << " " << name()
        << ": Sampling " << fieldMap_.size() << " fields" << endl;

    // Do analysis
    calcAddressing();

    return true;
}


Foam::wordList Foam::functionObjects::nearWallFields::fields() const
{
    wordList fields(fieldSet_.size());

    forAll(fieldSet_, fieldi)
    {
        fields[fieldi] = fieldSet_[fieldi].first();
    }

    return fields;
}


bool Foam::functionObjects::nearWallFields::execute()
{
    DebugInFunction << endl;

    if
    (
        fieldMap_.size()
     && vsf_.empty()
     && vvf_.empty()
     && vSpheretf_.empty()
     && vSymmtf_.empty()
     && vtf_.empty()
    )
    {
        Log << type() << " " << name()
            << ": Creating " << fieldMap_.size() << " fields" << endl;

        createFields(vsf_);
        createFields(vvf_);
        createFields(vSpheretf_);
        createFields(vSymmtf_);
        createFields(vtf_);

        Log << endl;
    }

    Log << type() << " " << name()
        << " write:" << nl
        << "    Sampling fields to " << time_.name()
        << endl;

    sampleFields(vsf_);
    sampleFields(vvf_);
    sampleFields(vSpheretf_);
    sampleFields(vSymmtf_);
    sampleFields(vtf_);

    return true;
}


bool Foam::functionObjects::nearWallFields::write()
{
    DebugInFunction << endl;

    Log << "    Writing sampled fields to " << time_.name()
        << endl;

    forAll(vsf_, i)
    {
        vsf_[i].write();
    }
    forAll(vvf_, i)
    {
        vvf_[i].write();
    }
    forAll(vSpheretf_, i)
    {
        vSpheretf_[i].write();
    }
    forAll(vSymmtf_, i)
    {
        vSymmtf_[i].write();
    }
    forAll(vtf_, i)
    {
        vtf_[i].write();
    }

    return true;
}


// ************************************************************************* //
