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

Description
    Best thing is probably to look at attachDetach which does almost exactly
    the same but for the geometric matching of points and face centres.

\*---------------------------------------------------------------------------*/

#include "perfectInterface.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "polyTopoChangeMap.H"
#include "matchPoints.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "polyRemoveFace.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(perfectInterface, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        perfectInterface,
        dictionary
    );
}


// Tolerance used as fraction of minimum edge length.
const Foam::scalar Foam::perfectInterface::tol_ = 1e-3;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointField Foam::perfectInterface::calcFaceCentres
(
    const indirectPrimitivePatch& pp
)
{
    const pointField& points = pp.points();

    pointField ctrs(pp.size());

    forAll(ctrs, patchFacei)
    {
        ctrs[patchFacei] = pp[patchFacei].centre(points);
    }

    return ctrs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::perfectInterface::perfectInterface
(
    const word& name,
    const label index,
    const polyTopoChanger& mme,
    const word& faceZoneName,
    const word& masterPatchName,
    const word& slavePatchName
)
:
    polyMeshModifier(name, index, mme, true),
    faceZoneID_(faceZoneName, mme.mesh().faceZones()),
    masterPatchID_(masterPatchName, mme.mesh().boundaryMesh()),
    slavePatchID_(slavePatchName, mme.mesh().boundaryMesh())
{}


// Construct from dictionary
Foam::perfectInterface::perfectInterface
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, readBool(dict.lookup("active"))),
    faceZoneID_
    (
        dict.lookup("faceZoneName"),
        mme.mesh().faceZones()
    ),
    masterPatchID_
    (
        dict.lookup("masterPatchName"),
        mme.mesh().boundaryMesh()
    ),
    slavePatchID_
    (
        dict.lookup("slavePatchName"),
        mme.mesh().boundaryMesh()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::perfectInterface::~perfectInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::perfectInterface::changeTopology() const
{
    // If modifier is inactive, skip change
    if (!active())
    {
        if (debug)
        {
            Pout<< "bool perfectInterface::changeTopology() const "
                << "for object " << name() << " : "
                << "Inactive" << endl;
        }

        return false;
    }
    else
    {
        // I want topo change every time step.
        return true;
    }
}


void Foam::perfectInterface::setRefinement
(
    const indirectPrimitivePatch& pp0,
    const indirectPrimitivePatch& pp1,
    polyTopoChange& ref
) const
{
    const polyMesh& mesh = topoChanger().mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Some aliases
    const edgeList& edges0 = pp0.edges();
    const pointField& pts0 = pp0.localPoints();
    const pointField& pts1 = pp1.localPoints();
    const labelList& meshPts0 = pp0.meshPoints();
    const labelList& meshPts1 = pp1.meshPoints();


    // Get local dimension as fraction of minimum edge length

    scalar minLen = great;

    forAll(edges0, edgeI)
    {
        minLen = min(minLen, edges0[edgeI].mag(pts0));
    }
    scalar typDim = tol_*minLen;

    if (debug)
    {
        Pout<< "typDim:" << typDim << " edges0:" << edges0.size()
            << " pts0:" << pts0.size() << " pts1:" << pts1.size()
            << " pp0:" << pp0.size() << " pp1:" << pp1.size() << endl;
    }


    // Determine pointMapping in mesh point labels. Uses geometric
    // comparison to find correspondence between patch points.

    labelList renumberPoints(mesh.points().size());
    forAll(renumberPoints, i)
    {
        renumberPoints[i] = i;
    }
    {
        labelList from1To0Points(pts1.size());

        bool matchOk = matchPoints
        (
            pts1,
            pts0,
            scalarField(pts1.size(), typDim),   // tolerance
            true,                               // verbose
            from1To0Points
        );

        if (!matchOk)
        {
            FatalErrorInFunction
                << "Points on patch sides do not match to within tolerance "
                << typDim << exit(FatalError);
        }

        forAll(pts1, i)
        {
            renumberPoints[meshPts1[i]] = meshPts0[from1To0Points[i]];
        }
    }



    // Calculate correspondence between patch faces

    labelList from0To1Faces(pp1.size());

    bool matchOk = matchPoints
    (
        calcFaceCentres(pp0),
        calcFaceCentres(pp1),
        scalarField(pp0.size(), typDim),    // tolerance
        true,                               // verbose
        from0To1Faces
    );

    if (!matchOk)
    {
        FatalErrorInFunction
            << "Face centres of patch sides do not match to within tolerance "
            << typDim << exit(FatalError);
    }



    // Now
    // - renumber faces using pts1 (except patch1 faces)
    // - remove patch1 faces. Remember cell label on owner side.
    // - modify patch0 faces to be internal.

    // 1. Get faces to be renumbered
    labelHashSet affectedFaces(2*pp1.size());
    forAll(meshPts1, i)
    {
        label meshPointi = meshPts1[i];

        if (meshPointi != renumberPoints[meshPointi])
        {
            const labelList& pFaces = mesh.pointFaces()[meshPointi];

            forAll(pFaces, pFacei)
            {
                affectedFaces.insert(pFaces[pFacei]);
            }
        }
    }
    forAll(pp1, i)
    {
        affectedFaces.erase(pp1.addressing()[i]);
    }
    // Remove patch0 from renumbered faces. Should not be necessary since
    // patch0 and 1 should not share any point (if created by mergeMeshing)
    // so affectedFaces should not contain any patch0 faces but you can
    // never be sure what the user is doing.
    forAll(pp0, i)
    {
        label facei = pp0.addressing()[i];

        if (affectedFaces.erase(facei))
        {
            WarningInFunction
                << "Found face " << facei << " vertices "
                << mesh.faces()[facei] << " whose points are"
                << " used both by master patch and slave patch" << endl;
        }
    }


    // 2. Renumber (non patch0/1) faces.
    forAllConstIter(labelHashSet, affectedFaces, iter)
    {
        const label facei = iter.key();
        const face& f = mesh.faces()[facei];

        face newFace(f.size());

        forAll(newFace, fp)
        {
            newFace[fp] = renumberPoints[f[fp]];
        }

        label nbr = -1;

        label patchi = -1;

        if (mesh.isInternalFace(facei))
        {
            nbr = mesh.faceNeighbour()[facei];
        }
        else
        {
            patchi = patches.whichPatch(facei);
        }

        label zoneID = mesh.faceZones().whichZone(facei);

        bool zoneFlip = false;

        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh.faceZones()[zoneID];

            zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
        }

        ref.setAction
        (
            polyModifyFace
            (
                newFace,                    // modified face
                facei,                      // label of face being modified
                mesh.faceOwner()[facei],    // owner
                nbr,                        // neighbour
                false,                      // face flip
                patchi,                     // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }


    // 3. Remove patch1 points
    forAll(meshPts1, i)
    {
        label meshPointi = meshPts1[i];

        if (meshPointi != renumberPoints[meshPointi])
        {
            ref.setAction(polyRemovePoint(meshPointi));
        }
    }


    // 4. Remove patch1 faces
    forAll(pp1, i)
    {
        label facei = pp1.addressing()[i];
        ref.setAction(polyRemoveFace(facei));
    }


    // 5. Modify patch0 faces for new points (not really necessary; see
    // comment above about patch1 and patch0 never sharing points) and
    // becoming internal.
    const boolList& mfFlip =
        mesh.faceZones()[faceZoneID_.index()].flipMap();

    forAll(pp0, i)
    {
        label facei = pp0.addressing()[i];

        const face& f = mesh.faces()[facei];

        face newFace(f.size());

        forAll(newFace, fp)
        {
            newFace[fp] = renumberPoints[f[fp]];
        }

        label own = mesh.faceOwner()[facei];

        label pp1Facei = pp1.addressing()[from0To1Faces[i]];

        label nbr = mesh.faceOwner()[pp1Facei];

        if (own < nbr)
        {
            ref.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    facei,                  // label of face being modified
                    own,                    // owner
                    nbr,                    // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    faceZoneID_.index(),    // zone for face
                    mfFlip[i]               // face flip in zone
                )
            );
        }
        else
        {
            ref.setAction
            (
                polyModifyFace
                (
                    newFace.reverseFace(),  // modified face
                    facei,                  // label of face being modified
                    nbr,                    // owner
                    own,                    // neighbour
                    true,                   // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    faceZoneID_.index(),    // zone for face
                    !mfFlip[i]              // face flip in zone
                )
            );
        }
    }
}


void Foam::perfectInterface::setRefinement(polyTopoChange& ref) const
{
    if (debug)
    {
        Pout<< "bool perfectInterface::setRefinement(polyTopoChange&) const : "
            << "for object " << name() << " : "
            << "masterPatchID_:" << masterPatchID_
            << " slavePatchID_:" << slavePatchID_
            << " faceZoneID_:" << faceZoneID_ << endl;
    }

    if
    (
        masterPatchID_.active()
     && slavePatchID_.active()
     && faceZoneID_.active()
    )
    {
        const polyMesh& mesh = topoChanger().mesh();

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const polyPatch& patch0 = patches[masterPatchID_.index()];
        const polyPatch& patch1 = patches[slavePatchID_.index()];


        labelList pp0Labels(identityMap(patch0.size())+patch0.start());
        indirectPrimitivePatch pp0
        (
            IndirectList<face>(mesh.faces(), pp0Labels),
            mesh.points()
        );

        labelList pp1Labels(identityMap(patch1.size())+patch1.start());
        indirectPrimitivePatch pp1
        (
            IndirectList<face>(mesh.faces(), pp1Labels),
            mesh.points()
        );

        setRefinement(pp0, pp1, ref);
    }
}


void Foam::perfectInterface::modifyMotionPoints(pointField& motionPoints) const
{
    // Update only my points. Nothing to be done here as points already
    // shared by now.
}


void Foam::perfectInterface::topoChange(const polyTopoChangeMap& map)
{
    // Mesh has changed topologically.  Update local topological data
    const polyMesh& mesh = topoChanger().mesh();

    faceZoneID_.update(mesh.faceZones());
    masterPatchID_.update(mesh.boundaryMesh());
    slavePatchID_.update(mesh.boundaryMesh());
}


void Foam::perfectInterface::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name()<< nl
        << faceZoneID_.name() << nl
        << masterPatchID_.name() << nl
        << slavePatchID_.name() << endl;
}


void Foam::perfectInterface::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl

        << "    type " << type()
        << token::END_STATEMENT << nl

        << "    active " << active()
        << token::END_STATEMENT << nl

        << "    faceZoneName " << faceZoneID_.name()
        << token::END_STATEMENT << nl

        << "    masterPatchName " << masterPatchID_.name()
        << token::END_STATEMENT << nl

        << "    slavePatchName " << slavePatchID_.name()
        << token::END_STATEMENT << nl

        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
