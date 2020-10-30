/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "extrudePatchMesh.H"

#include "createShellMesh.H"
#include "polyTopoChange.H"
#include "wallPolyPatch.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(extrudePatchMesh, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrudePatchMesh::extrudePatchMesh
(
    const fvMesh& mesh,
    const fvPatch& patch,
    const dictionary& dict,
    const word regionName,
    const List<polyPatch*>& regionPatches
)
:
    fvMesh
    (
        IOobject
        (
            regionName,
            mesh.facesInstance(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true
        ),
        pointField(),
        faceList(),
        labelList(),
        labelList(),
        false
    ),
    extrudedPatch_(patch.patch()),
    dict_(dict)
{
    extrudeMesh(regionPatches);
}


extrudePatchMesh::extrudePatchMesh
(
    const fvMesh& mesh,
    const fvPatch& patch,
    const dictionary& dict,
    const word regionName
)
:
    fvMesh
    (
        IOobject
        (
            regionName,
            mesh.facesInstance(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true
        ),
        pointField(),
        faceList(),
        labelList(),
        labelList(),
        false
    ),
    extrudedPatch_(patch.patch()),
    dict_(dict)
{

    List<polyPatch*> regionPatches(3);
    List<word> patchNames(regionPatches.size());
    List<word> patchTypes(regionPatches.size());
    PtrList<dictionary> dicts(regionPatches.size());

    forAll(dicts, patchi)
    {
        if (!dicts.set(patchi))
        {
            dicts.set(patchi, new dictionary());
        }
    }

    dicts[bottomPatchID] = dict_.subDict("bottomCoeffs");
    dicts[sidePatchID] = dict_.subDict("sideCoeffs");
    dicts[topPatchID] = dict_.subDict("topCoeffs");

    forAll(dicts, patchi)
    {
        dicts[patchi].lookup("name") >> patchNames[patchi];
        dicts[patchi].lookup("type") >> patchTypes[patchi];
    }

    forAll(regionPatches, patchi)
    {
        dictionary&  patchDict = dicts[patchi];
        patchDict.set("nFaces", 0);
        patchDict.set("startFace", 0);

        regionPatches[patchi] = polyPatch::New
            (
                patchNames[patchi],
                patchDict,
                patchi,
                mesh.boundaryMesh()
            ).ptr();

    }

    extrudeMesh(regionPatches);

}


void extrudePatchMesh::extrudeMesh(const List<polyPatch*>& regionPatches)
{
    if (this->boundaryMesh().size() == 0)
    {
        bool columnCells = readBool(dict_.lookup("columnCells"));

        PackedBoolList nonManifoldEdge(extrudedPatch_.nEdges());
        for (label edgeI = 0; edgeI < extrudedPatch_.nInternalEdges(); edgeI++)
        {
            if (columnCells)
            {
                nonManifoldEdge[edgeI] = true;
            }
        }

        autoPtr<extrudeModel> model_(extrudeModel::New(dict_));

        faceList pointGlobalRegions;
        faceList pointLocalRegions;
        labelList localToGlobalRegion;

        const primitiveFacePatch pp
        (
            extrudedPatch_, extrudedPatch_.points()
        );

        createShellMesh::calcPointRegions
        (
            this->globalData(),
            pp,
            nonManifoldEdge,
            false,

            pointGlobalRegions,
            pointLocalRegions,
            localToGlobalRegion
        );


        // Per local region an originating point
        labelList localRegionPoints(localToGlobalRegion.size());
        forAll(pointLocalRegions, facei)
        {
            const face& f = extrudedPatch_.localFaces()[facei];
            const face& pRegions = pointLocalRegions[facei];
            forAll(pRegions, fp)
            {
                localRegionPoints[pRegions[fp]] = f[fp];
            }
        }

        // Calculate region normals by reducing local region normals
        pointField localRegionNormals(localToGlobalRegion.size());
        {
            pointField localSum(localToGlobalRegion.size(), Zero);

            forAll(pointLocalRegions, facei)
            {
                const face& pRegions = pointLocalRegions[facei];
                forAll(pRegions, fp)
                {
                    label localRegionI = pRegions[fp];
                    localSum[localRegionI] +=
                        extrudedPatch_.faceNormals()[facei];
                }
            }

            Map<point> globalSum(2*localToGlobalRegion.size());

            forAll(localSum, localRegionI)
            {
                label globalRegionI = localToGlobalRegion[localRegionI];
                globalSum.insert(globalRegionI, localSum[localRegionI]);
            }

            // Reduce
            Pstream::mapCombineGather(globalSum, plusEqOp<point>());
            Pstream::mapCombineScatter(globalSum);

            forAll(localToGlobalRegion, localRegionI)
            {
                label globalRegionI = localToGlobalRegion[localRegionI];
                localRegionNormals[localRegionI] = globalSum[globalRegionI];
            }
            localRegionNormals /= mag(localRegionNormals);
        }


        // Per local region an extrusion direction
        vectorField firstDisp(localToGlobalRegion.size());
        forAll(firstDisp, regionI)
        {
            // const point& regionPt = regionCentres[regionI];
            const point& regionPt = extrudedPatch_.points()
            [
                extrudedPatch_.meshPoints()
                [
                    localRegionPoints[regionI]
                ]
            ];
            const vector& n = localRegionNormals[regionI];
            firstDisp[regionI] = model_()(regionPt, n, 1) - regionPt;
        }


        // Extrude engine
        createShellMesh extruder
        (
            pp,
            pointLocalRegions,
            localRegionPoints
        );
/*
        List<polyPatch*> regionPatches(3);
        List<word> patchNames(regionPatches.size());
        List<word> patchTypes(regionPatches.size());
        PtrList<dictionary> dicts(regionPatches.size());

        forAll(dicts, patchi)
        {
            if (!dicts.set(patchi))
            {
                dicts.set(patchi, new dictionary());
            }
        }

        dicts[bottomPatchID] = dict_.subDict("bottomCoeffs");
        dicts[sidePatchID] = dict_.subDict("sideCoeffs");
        dicts[topPatchID] = dict_.subDict("topCoeffs");

        forAll(dicts, patchi)
        {
            dicts[patchi].lookup("name") >> patchNames[patchi];
            dicts[patchi].lookup("type") >> patchTypes[patchi];
        }

        forAll(regionPatches, patchi)
        {
            dictionary&  patchDict = dicts[patchi];
            patchDict.set("nFaces", 0);
            patchDict.set("startFace", 0);

            regionPatches[patchi] = polyPatch::New
                (
                    patchNames[patchi],
                    patchDict,
                    patchi,
                    mesh.boundaryMesh()
                ).ptr();

        }
*/
        this->clearOut();
        this->removeFvBoundary();
        this->addFvPatches(regionPatches, true);


        // At this point we have a valid mesh with 3 patches and zero cells.
        // Determine:
        // - per face the top and bottom patch (topPatchID, bottomPatchID)
        // - per edge, per face on edge the side patch (edgePatches)
        labelListList edgePatches(extrudedPatch_.nEdges());
        forAll(edgePatches, edgeI)
        {
            const labelList& eFaces = extrudedPatch_.edgeFaces()[edgeI];

            if (eFaces.size() != 2 || nonManifoldEdge[edgeI])
            {
                edgePatches[edgeI].setSize(eFaces.size(), sidePatchID);
            }
        }

        polyTopoChange meshMod(regionPatches.size());

        extruder.setRefinement
        (
            firstDisp,                              // first displacement
            model_().expansionRatio(),
            model_().nLayers(),                       // nLayers
            labelList(extrudedPatch_.size(), topPatchID),
            labelList(extrudedPatch_.size(), bottomPatchID),
            edgePatches,
            meshMod
        );

        autoPtr<mapPolyMesh> map = meshMod.changeMesh
        (
            *this,              // mesh to change
            false               // inflate
        );

        // Update numbering on extruder.
        extruder.updateMesh(map);

        this->setInstance(this->thisDb().time().constant());
        this->write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

extrudePatchMesh::~extrudePatchMesh()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
