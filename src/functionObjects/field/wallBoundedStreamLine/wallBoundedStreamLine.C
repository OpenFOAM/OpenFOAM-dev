/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "wallBoundedStreamLine.H"
#include "fvMesh.H"
#include "wallBoundedStreamLineParticleCloud.H"
#include "ReadFields.H"
#include "meshSearch.H"
#include "sampledSet.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "interpolationCellPoint.H"
#include "PatchTools.H"
#include "meshSearchMeshObject.H"
#include "faceSet.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallBoundedStreamLine, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        wallBoundedStreamLine,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::functionObjects::wallBoundedStreamLine::wallPatch() const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nFaces = 0;

    forAll(patches, patchi)
    {
        //if (!polyPatch::constraintType(patches[patchi].type()))
        if (isA<wallPolyPatch>(patches[patchi]))
        {
            nFaces += patches[patchi].size();
        }
    }

    labelList addressing(nFaces);

    nFaces = 0;

    forAll(patches, patchi)
    {
        //if (!polyPatch::constraintType(patches[patchi].type()))
        if (isA<wallPolyPatch>(patches[patchi]))
        {
            const polyPatch& pp = patches[patchi];

            forAll(pp, i)
            {
                addressing[nFaces++] = pp.start()+i;
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                addressing
            ),
            mesh.points()
        )
    );
}


Foam::tetIndices Foam::functionObjects::wallBoundedStreamLine::findNearestTet
(
    const PackedBoolList& isWallPatch,
    const point& seedPt,
    const label celli
) const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    const cell& cFaces = mesh.cells()[celli];

    label minFacei = -1;
    label minTetPtI = -1;
    scalar minDistSqr = sqr(GREAT);

    forAll(cFaces, cFacei)
    {
        label facei = cFaces[cFacei];

        if (isWallPatch[facei])
        {
            const face& f = mesh.faces()[facei];
            const label fp0 = mesh.tetBasePtIs()[facei];
            const point& basePoint = mesh.points()[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = mesh.points()[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = mesh.points()[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);

                scalar d2 = magSqr(tri.centre() - seedPt);
                if (d2 < minDistSqr)
                {
                    minDistSqr = d2;
                    minFacei = facei;
                    minTetPtI = i-1;
                }
                fp = nextFp;
            }
        }
    }

    // Put particle in tet
    return tetIndices
    (
        celli,
        minFacei,
        minTetPtI,
        mesh
    );
}


void Foam::functionObjects::wallBoundedStreamLine::track()
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);


    // Determine the 'wall' patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // These are the faces that need to be followed

    autoPtr<indirectPrimitivePatch> boundaryPatch(wallPatch());
    PackedBoolList isWallPatch(mesh.nFaces());
    forAll(boundaryPatch().addressing(), i)
    {
        isWallPatch[boundaryPatch().addressing()[i]] = 1;
    }



    // Find nearest wall particle for the seedPoints
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IDLList<wallBoundedStreamLineParticle> initialParticles;
    wallBoundedStreamLineParticleCloud particles
    (
        mesh,
        cloudName_,
        initialParticles
    );

    {
        // Get the seed points
        // ~~~~~~~~~~~~~~~~~~~

        const sampledSet& seedPoints = sampledSetPtr_();


        forAll(seedPoints, i)
        {
            const point& seedPt = seedPoints[i];
            label celli = seedPoints.cells()[i];

            tetIndices ids(findNearestTet(isWallPatch, seedPt, celli));

            if (ids.face() != -1 && isWallPatch[ids.face()])
            {
                //Pout<< "Seeding particle :" << nl
                //    << "     seedPt:" << seedPt << nl
                //    << "     face  :" << ids.face() << nl
                //    << "     at    :" << mesh.faceCentres()[ids.face()] << nl
                //    << "     cell  :" << mesh.cellCentres()[ids.cell()] << nl
                //    << endl;

                particles.addParticle
                (
                    new wallBoundedStreamLineParticle
                    (
                        mesh,
                        ids.faceTri(mesh).centre(),
                        ids.cell(),
                        ids.face(),     // tetFace
                        ids.tetPt(),
                        -1,             // not on a mesh edge
                        -1,             // not on a diagonal edge
                        lifeTime_       // lifetime
                    )
                );
            }
            else
            {
                Pout<< type() << " : ignoring seed " << seedPt
                    << " since not in wall cell." << endl;
            }
        }
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    Info<< type() << " : seeded " << nSeeds << " particles." << endl;



    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolation<scalar>> vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolation<vector>> vvInterp;

    label UIndex = -1;

    label nScalar = 0;
    label nVector = 0;

    forAll(fields_, i)
    {
        if (mesh.foundObject<volScalarField>(fields_[i]))
        {
            nScalar++;
        }
        else if (mesh.foundObject<volVectorField>(fields_[i]))
        {
            nVector++;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find field " << fields_[i] << endl
                << "Valid scalar fields are:"
                << mesh.names(volScalarField::typeName) << endl
                << "Valid vector fields are:"
                << mesh.names(volVectorField::typeName)
                << exit(FatalError);
        }
    }

    vsInterp.setSize(nScalar);
    nScalar = 0;
    vvInterp.setSize(nVector);
    nVector = 0;

    forAll(fields_, i)
    {
        if (mesh.foundObject<volScalarField>(fields_[i]))
        {
            const volScalarField& f = mesh.lookupObject<volScalarField>
            (
                fields_[i]
            );
            vsInterp.set
            (
                nScalar++,
                interpolation<scalar>::New
                (
                    interpolationScheme_,
                    f
                )
            );
        }
        else if (mesh.foundObject<volVectorField>(fields_[i]))
        {
            const volVectorField& f = mesh.lookupObject<volVectorField>
            (
                fields_[i]
            );

            if (f.name() == UName_)
            {
                UIndex = nVector;
            }

            vvInterp.set
            (
                nVector++,
                interpolation<vector>::New
                (
                    interpolationScheme_,
                    f
                )
            );
        }
    }

    // Store the names
    scalarNames_.setSize(vsInterp.size());
    forAll(vsInterp, i)
    {
        scalarNames_[i] = vsInterp[i].psi().name();
    }
    vectorNames_.setSize(vvInterp.size());
    forAll(vvInterp, i)
    {
        vectorNames_[i] = vvInterp[i].psi().name();
    }

    // Check that we know the index of U in the interpolators.

    if (UIndex == -1)
    {
        FatalErrorInFunction
            << "Cannot find field to move particles with : " << UName_
            << endl
            << "This field has to be present in the sampled fields "
            << fields_
            << " and in the objectRegistry." << endl
            << exit(FatalError);
    }

    // Sampled data
    // ~~~~~~~~~~~~

    // Size to maximum expected sizes.
    allTracks_.clear();
    allTracks_.setCapacity(nSeeds);
    allScalars_.setSize(vsInterp.size());
    forAll(allScalars_, i)
    {
        allScalars_[i].clear();
        allScalars_[i].setCapacity(nSeeds);
    }
    allVectors_.setSize(vvInterp.size());
    forAll(allVectors_, i)
    {
        allVectors_[i].clear();
        allVectors_[i].setCapacity(nSeeds);
    }


    // additional particle info
    wallBoundedStreamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction?
        trackLength_,   // fixed track length
        isWallPatch,    // which faces are to follow

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoundedStreamLine::wallBoundedStreamLine
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    dict_(dict)
{
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoundedStreamLine::~wallBoundedStreamLine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallBoundedStreamLine::read(const dictionary& dict)
{
    //dict_ = dict;
    dict.lookup("fields") >> fields_;
    if (dict.found("U"))
    {
        dict.lookup("U") >> UName_;
    }
    else
    {
        UName_ = "U";
        if (dict.found("U"))
        {
            IOWarningInFunction
            (
                dict
            )   << "Using deprecated entry \"U\"."
                << " Please use \"UName\" instead."
                << endl;
            dict.lookup("U") >> UName_;
        }
    }

    if (findIndex(fields_, UName_) == -1)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Velocity field for tracking " << UName_
            << " should be present in the list of fields " << fields_
            << exit(FatalIOError);
    }


    dict.lookup("trackForward") >> trackForward_;
    dict.lookup("lifeTime") >> lifeTime_;
    if (lifeTime_ < 1)
    {
        FatalErrorInFunction
            << "Illegal value " << lifeTime_ << " for lifeTime"
            << exit(FatalError);
    }
    trackLength_ = VGREAT;
    if (dict.found("trackLength"))
    {
        dict.lookup("trackLength") >> trackLength_;

        Info<< type() << " : fixed track length specified : "
            << trackLength_ << nl << endl;
    }


    interpolationScheme_ = dict.lookupOrDefault
    (
        "interpolationScheme",
        interpolationCellPoint<scalar>::typeName
    );

    cloudName_ = dict.lookupOrDefault<word>
    (
        "cloudName",
        "wallBoundedStreamLine"
    );
    dict.lookup("seedSampleSet") >> seedSet_;

    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    const dictionary& coeffsDict = dict.subDict(seedSet_ + "Coeffs");

    sampledSetPtr_ = sampledSet::New
    (
        seedSet_,
        mesh,
        meshSearchMeshObject::New(mesh),
        coeffsDict
    );
    coeffsDict.lookup("axis") >> sampledSetAxis_;

    scalarFormatterPtr_ = writer<scalar>::New(dict.lookup("setFormat"));
    vectorFormatterPtr_ = writer<vector>::New(dict.lookup("setFormat"));


    // Make sure that the mesh is trackable
    if (debug)
    {
        // 1. positive volume decomposition tets
        faceSet faces(mesh, "lowQualityTetFaces", mesh.nFaces()/100+1);
        if
        (
            polyMeshTetDecomposition::checkFaceTets
            (
                mesh,
                polyMeshTetDecomposition::minTetQuality,
                true,
                &faces
            )
        )
        {
            label nFaces = returnReduce(faces.size(), sumOp<label>());

            WarningInFunction
                << "Found " << nFaces
                <<" faces with low quality or negative volume "
                << "decomposition tets. Writing to faceSet " << faces.name()
                << endl;
        }

        // 2. all edges on a cell having two faces
        EdgeMap<label> numFacesPerEdge;
        forAll(mesh.cells(), celli)
        {
            const cell& cFaces = mesh.cells()[celli];

            numFacesPerEdge.clear();

            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];
                const face& f = mesh.faces()[facei];
                forAll(f, fp)
                {
                    const edge e(f[fp], f.nextLabel(fp));
                    EdgeMap<label>::iterator eFnd =
                        numFacesPerEdge.find(e);
                    if (eFnd != numFacesPerEdge.end())
                    {
                        eFnd()++;
                    }
                    else
                    {
                        numFacesPerEdge.insert(e, 1);
                    }
                }
            }

            forAllConstIter(EdgeMap<label>, numFacesPerEdge, iter)
            {
                if (iter() != 2)
                {
                    FatalErrorInFunction
                        << "problem cell:" << celli
                        << abort(FatalError);
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::wallBoundedStreamLine::execute()
{
    return true;
}


bool Foam::functionObjects::wallBoundedStreamLine::write()
{
    const Time& runTime = obr_.time();
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);


    // Do all injection and tracking
    track();


    if (Pstream::parRun())
    {
        // Append slave tracks to master ones
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        globalIndex globalTrackIDs(allTracks_.size());

        // Construct a distribution map to pull all to the master.
        labelListList sendMap(Pstream::nProcs());
        labelListList recvMap(Pstream::nProcs());

        if (Pstream::master())
        {
            // Master: receive all. My own first, then consecutive
            // processors.
            label trackI = 0;

            forAll(recvMap, proci)
            {
                labelList& fromProc = recvMap[proci];
                fromProc.setSize(globalTrackIDs.localSize(proci));
                forAll(fromProc, i)
                {
                    fromProc[i] = trackI++;
                }
            }
        }

        labelList& toMaster = sendMap[0];
        toMaster.setSize(globalTrackIDs.localSize());
        forAll(toMaster, i)
        {
            toMaster[i] = i;
        }

        const mapDistribute distMap
        (
            globalTrackIDs.size(),
            sendMap.xfer(),
            recvMap.xfer()
        );

        // Distribute the track positions. Note: use scheduled comms
        // to prevent buffering.
        allTracks_.shrink();
        mapDistributeBase::distribute
        (
            Pstream::scheduled,
            distMap.schedule(),
            distMap.constructSize(),
            distMap.subMap(),
            false,
            distMap.constructMap(),
            false,
            allTracks_,
            flipOp()
        );

        // Distribute the scalars
        forAll(allScalars_, scalarI)
        {
            allScalars_[scalarI].shrink();
            mapDistributeBase::distribute
            (
                Pstream::scheduled,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                false,
                distMap.constructMap(),
                false,
                allScalars_[scalarI],
                flipOp()
            );
            allScalars_[scalarI].setCapacity(allScalars_[scalarI].size());
        }
        // Distribute the vectors
        forAll(allVectors_, vectorI)
        {
            allVectors_[vectorI].shrink();
            mapDistributeBase::distribute
            (
                Pstream::scheduled,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),
                false,
                distMap.constructMap(),
                false,
                allVectors_[vectorI],
                flipOp()
            );
            allVectors_[vectorI].setCapacity(allVectors_[vectorI].size());
        }
    }


    label n = 0;
    forAll(allTracks_, trackI)
    {
        n += allTracks_[trackI].size();
    }

    Info<< "    Tracks:" << allTracks_.size() << nl
        << "    Total samples:" << n << endl;


    // Massage into form suitable for writers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (Pstream::master() && allTracks_.size())
    {
        // Make output directory

        fileName vtkPath
        (
            Pstream::parRun()
          ? runTime.path()/".."/"postProcessing"/"sets"/name()
          : runTime.path()/"postProcessing"/"sets"/name()
        );
        if (mesh.name() != fvMesh::defaultRegion)
        {
            vtkPath = vtkPath/mesh.name();
        }
        vtkPath = vtkPath/mesh.time().timeName();

        mkDir(vtkPath);

        // Convert track positions

        PtrList<coordSet> tracks(allTracks_.size());
        forAll(allTracks_, trackI)
        {
            tracks.set
            (
                trackI,
                new coordSet
                (
                    "track" + Foam::name(trackI),
                    sampledSetAxis_                 //"xyz"
                )
            );
            tracks[trackI].transfer(allTracks_[trackI]);
        }

        // Convert scalar values

        if (allScalars_.size() > 0)
        {
            List<List<scalarField>> scalarValues(allScalars_.size());

            forAll(allScalars_, scalarI)
            {
                DynamicList<scalarList>& allTrackVals =
                    allScalars_[scalarI];
                scalarValues[scalarI].setSize(allTrackVals.size());

                forAll(allTrackVals, trackI)
                {
                    scalarList& trackVals = allTrackVals[trackI];
                    scalarValues[scalarI][trackI].transfer(trackVals);
                }
            }

            fileName vtkFile
            (
                vtkPath
              / scalarFormatterPtr_().getFileName
                (
                    tracks[0],
                    scalarNames_
                )
            );

            Info<< "    Writing data to " << vtkFile.path() << endl;

            scalarFormatterPtr_().write
            (
                true,           // writeTracks
                tracks,
                scalarNames_,
                scalarValues,
                OFstream(vtkFile)()
            );
        }

        // Convert vector values

        if (allVectors_.size() > 0)
        {
            List<List<vectorField>> vectorValues(allVectors_.size());

            forAll(allVectors_, vectorI)
            {
                DynamicList<vectorList>& allTrackVals =
                    allVectors_[vectorI];
                vectorValues[vectorI].setSize(allTrackVals.size());

                forAll(allTrackVals, trackI)
                {
                    vectorList& trackVals = allTrackVals[trackI];
                    vectorValues[vectorI][trackI].transfer(trackVals);
                }
            }

            fileName vtkFile
            (
                vtkPath
              / vectorFormatterPtr_().getFileName
                (
                    tracks[0],
                    vectorNames_
                )
            );

            vectorFormatterPtr_().write
            (
                true,           // writeTracks
                tracks,
                vectorNames_,
                vectorValues,
                OFstream(vtkFile)()
            );
        }
    }

    return true;
}


void Foam::functionObjects::wallBoundedStreamLine::updateMesh
(
    const mapPolyMesh& mpm
)
{
    const fvMesh& mesh_ = dynamic_cast<const fvMesh&>(obr_);

    if (&mpm.mesh() == &mesh_)
    {
        read(dict_);
    }
}


void Foam::functionObjects::wallBoundedStreamLine::movePoints
(
    const polyMesh& mesh
)
{
    const fvMesh& mesh_ = dynamic_cast<const fvMesh&>(obr_);

    if (&mesh == &mesh_)
    {
        // Moving mesh affects the search tree
        read(dict_);
    }
}


// ************************************************************************* //
