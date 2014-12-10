/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "mappedPatchBase.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearWallFields, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nearWallFields::calcAddressing()
{
    const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Count number of faces
    label nPatchFaces = 0;
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        nPatchFaces += mesh.boundary()[patchI].size();
    }

    // Global indexing
    globalIndex globalWalls(nPatchFaces);

    if (debug)
    {
        Info<< "nearWallFields::calcAddressing() :"
            << " nPatchFaces:" << globalWalls.size() << endl;
    }

    // Construct cloud
    Cloud<findCellParticle> cloud(mesh, IDLList<findCellParticle>());

    // Add particles to track to sample locations
    nPatchFaces = 0;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const fvPatch& patch = mesh.boundary()[patchI];

        vectorField nf(patch.nf());
        vectorField faceCellCentres(patch.patch().faceCellCentres());

        forAll(patch, patchFaceI)
        {
            label meshFaceI = patch.start()+patchFaceI;

            // Find starting point on face (since faceCentre might not
            // be on face-diagonal decomposition)
            pointIndexHit startInfo
            (
                mappedPatchBase::facePoint
                (
                    mesh,
                    meshFaceI,
                    polyMesh::FACEDIAGTETS
                )
            );


            point start;
            if (startInfo.hit())
            {
                start = startInfo.hitPoint();
            }
            else
            {
                // Fallback: start tracking from neighbouring cell centre
                start = faceCellCentres[patchFaceI];
            }

            const point end = start-distance_*nf[patchFaceI];

            // Find tet for starting location
            label cellI = -1;
            label tetFaceI = -1;
            label tetPtI = -1;
            mesh.findCellFacePt(start, cellI, tetFaceI, tetPtI);

            // Add to cloud. Add originating face as passive data
            cloud.addParticle
            (
                new findCellParticle
                (
                    mesh,
                    start,
                    cellI,
                    tetFaceI,
                    tetPtI,
                    end,
                    globalWalls.toGlobal(nPatchFaces)    // passive data
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
            mesh.time().path()
           /"wantedTracks_" + mesh.time().timeName() + ".obj"
        );
        Info<< "nearWallFields::calcAddressing() :"
            << "Dumping tracks to " << str.name() << endl;

        forAllConstIter(Cloud<findCellParticle>, cloud, iter)
        {
            const findCellParticle& tp = iter();
            str.write(linePointRef(tp.position(), tp.end()));
        }
    }



    // Per cell: empty or global wall index and end location
    cellToWalls_.setSize(mesh.nCells());
    cellToSamples_.setSize(mesh.nCells());

    // Database to pass into findCellParticle::move
    findCellParticle::trackingData td(cloud, cellToWalls_, cellToSamples_);

    // Track all particles to their end position.
    scalar maxTrackLen = 2.0*mesh.bounds().mag();


    //Debug: collect start points
    pointField start;
    if (debug)
    {
        start.setSize(nPatchFaces);
        nPatchFaces = 0;
        forAllConstIter(Cloud<findCellParticle>, cloud, iter)
        {
            const findCellParticle& tp = iter();
            start[nPatchFaces++] = tp.position();
        }
    }


    cloud.move(td, maxTrackLen);


    // Rework cell-to-globalpatchface into a map
    List<Map<label> > compactMap;
    getPatchDataMapPtr_.reset
    (
        new mapDistribute
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
                mesh.time().path()
               /"obtainedTracks_" + mesh.time().timeName() + ".obj"
            );
            Info<< "nearWallFields::calcAddressing() :"
                << "Dumping obtained to " << str.name() << endl;

            forAll(cellToWalls_, cellI)
            {
                const List<point>& ends = cellToSamples_[cellI];
                const labelList& cData = cellToWalls_[cellI];
                forAll(cData, i)
                {
                    str.write(linePointRef(ends[i], start[cData[i]]));
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearWallFields::nearWallFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    fieldSet_()
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "nearWallFields::nearWallFields"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_
            << endl;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nearWallFields::~nearWallFields()
{
    if (debug)
    {
        Info<< "nearWallFields::~nearWallFields()" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nearWallFields::read(const dictionary& dict)
{
    if (debug)
    {
        Info<< "nearWallFields::read(const dictionary&)" << endl;
    }

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        dict.lookup("fields") >> fieldSet_;
        patchSet_ =
            mesh.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));
        distance_ = readScalar(dict.lookup("distance"));


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
        forAll(fieldSet_, setI)
        {
            const word& fldName = fieldSet_[setI].first();
            const word& sampleFldName = fieldSet_[setI].second();

            fieldMap_.insert(fldName, sampleFldName);
            reverseFieldMap_.insert(sampleFldName, fldName);
        }

        Info<< type() << " " << name_ << ": Sampling " << fieldMap_.size()
            << " fields" << endl;

        // Do analysis
        calcAddressing();
    }
}


void Foam::nearWallFields::execute()
{
    if (debug)
    {
        Info<< "nearWallFields:execute()" << endl;
    }


    if (active_)
    {
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
            Info<< type() << " " << name_ << ": Creating " << fieldMap_.size()
                << " fields" << endl;

            createFields(vsf_);
            createFields(vvf_);
            createFields(vSpheretf_);
            createFields(vSymmtf_);
            createFields(vtf_);

            Info<< endl;
        }

        Info<< type() << " " << name_ << " output:" << nl;

        Info<< "    Sampling fields to " << obr_.time().timeName()
            << endl;

        sampleFields(vsf_);
        sampleFields(vvf_);
        sampleFields(vSpheretf_);
        sampleFields(vSymmtf_);
        sampleFields(vtf_);
    }
}


void Foam::nearWallFields::end()
{
    if (debug)
    {
        Info<< "nearWallFields:end()" << endl;
    }

    if (active_)
    {
        execute();
    }
}


void Foam::nearWallFields::timeSet()
{
    // Do nothing
}


void Foam::nearWallFields::write()
{
    if (debug)
    {
        Info<< "nearWallFields:write()" << endl;
    }

    if (active_)
    {
        Info<< "    Writing sampled fields to " << obr_.time().timeName()
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

        Info<< endl;
    }
}


// ************************************************************************* //
