/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

Application
    testPatchTools

Description
    Test app for PatchTools functionality

\*---------------------------------------------------------------------------*/

#include "PatchTools.H"
#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
#include "OBJstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//template<class PatchType>
//Foam::tmp<Foam::pointField>
//areaPointNormals
//(
//    const polyMesh& mesh,
//    const PatchType& p,
//    const labelList& meshFaces
//)
//{
//    // Assume patch is smaller than the globalData().coupledPatch() (?) so
//    // loop over patch meshPoints.
//
//    const labelList& meshPoints = p.meshPoints();
//
//    const globalMeshData& globalData = mesh.globalData();
//    const indirectPrimitivePatch& coupledPatch = globalData.coupledPatch();
//    const Map<label>& coupledPatchMP = coupledPatch.meshPointMap();
//    const mapDistribute& map = globalData.globalPointSlavesMap();
//    const globalIndexAndTransform& transforms =
//        globalData.globalTransforms();
//
//
//    // 1. Start off with local (area-weighted) normals
//    //    (note:without calculating pointNormals
//    //     to avoid them being stored)
//
//    tmp<pointField> textrudeN(new pointField(p.nPoints(), Zero));
//    pointField& extrudeN = textrudeN();
//    {
//        const faceList& localFaces = p.localFaces();
//        const vectorField& faceAreas = mesh.faceAreas();
//
//        forAll(localFaces, facei)
//        {
//            const face& f = localFaces[facei];
//            const vector& n = faceAreas[meshFaces[facei]];
//            forAll(f, fp)
//            {
//                extrudeN[f[fp]] += n;
//            }
//        }
//    }
//
//
//    // Collect local pointFaces
//    List<List<point>> pointFaceNormals(map.constructSize());
//    {
//        const vectorField& faceAreas = mesh.faceAreas();
//
//        forAll(meshPoints, patchPointi)
//        {
//            label meshPointi = meshPoints[patchPointi];
//            Map<label>::const_iterator fnd = coupledPatchMP.find(meshPointi);
//            if (fnd != coupledPatchMP.end())
//            {
//                label coupledPointi = fnd();
//
//                List<point>& pNormals = pointFaceNormals[coupledPointi];
//                const labelList& pFaces = p.pointFaces()[patchPointi];
//                pNormals.setSize(pFaces.size());
//                forAll(pFaces, i)
//                {
//                    pNormals[i] =  faceAreas[meshFaces[pFaces[i]]];
//                }
//            }
//        }
//    }
//
//    // Pull remote data into local slots
//    map.distribute
//    (
//        transforms,
//        pointFaceNormals,
//        listTransform()
//    );
//
//
//    // Combine normals
//    const labelListList& slaves = globalData.globalPointSlaves();
//    const labelListList& transformedSlaves =
//        globalData.globalPointTransformedSlaves();
//
//
//    pointField coupledPointNormals(map.constructSize(), Zero);
//
//    forAll(meshPoints, patchPointi)
//    {
//        label meshPointi = meshPoints[patchPointi];
//        Map<label>::const_iterator fnd = coupledPatchMP.find(meshPointi);
//        if (fnd != coupledPatchMP.end())
//        {
//            label coupledPointi = fnd();
//            const labelList& slaveSlots = slaves[coupledPointi];
//            const labelList& transformedSlaveSlots =
//                transformedSlaves[coupledPointi];
//
//            label nFaces = slaveSlots.size()+transformedSlaveSlots.size();
//            if (nFaces > 0)
//            {
//                // Combine
//                point& n = coupledPointNormals[coupledPointi];
//
//                n += sum(pointFaceNormals[coupledPointi]);
//
//                forAll(slaveSlots, i)
//                {
//                    n += sum(pointFaceNormals[slaveSlots[i]]);
//                }
//                forAll(transformedSlaveSlots, i)
//                {
//                    n += sum(pointFaceNormals[transformedSlaveSlots[i]]);
//                }
//
//                // Put back into slave slots
//                forAll(slaveSlots, i)
//                {
//                    coupledPointNormals[slaveSlots[i]] = n;
//                }
//                forAll(transformedSlaveSlots, i)
//                {
//                    coupledPointNormals[transformedSlaveSlots[i]] = n;
//                }
//            }
//        }
//    }
//
//
//    // Send back
//    map.reverseDistribute
//    (
//        transforms,
//        coupledPointNormals.size(),
//        coupledPointNormals,
//        mapDistribute::transform()
//    );
//
//
//    // Override patch normals
//    forAll(meshPoints, patchPointi)
//    {
//        label meshPointi = meshPoints[patchPointi];
//        Map<label>::const_iterator fnd = coupledPatchMP.find(meshPointi);
//        if (fnd != coupledPatchMP.end())
//        {
//            label coupledPointi = fnd();
//            extrudeN[patchPointi] = coupledPointNormals[coupledPointi];
//        }
//    }
//
//    extrudeN /= mag(extrudeN)+vSmall;
//
//    return textrudeN;
//}



// Main program:

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    argList::validArgs.append("patch");
    #include "setRootCase.H"
    #include "createTime.H"

    #include "createMesh.H"

    const word patchName = args[1];
    label patchi = mesh.boundaryMesh().findPatchID(patchName);
    const polyPatch& pp = mesh.boundaryMesh()[patchi];

    const indirectPrimitivePatch& cpp = mesh.globalData().coupledPatch();

    {
        OBJstream str(runTime.path()/"edgePatchNormals.obj");

        labelList patchEdges;
        labelList coupledEdges;
        PackedBoolList sameEdgeOrientation;
        PatchTools::matchEdges
        (
            pp,
            cpp,
            patchEdges,
            coupledEdges,
            sameEdgeOrientation
        );

        const pointField en
        (
            PatchTools::edgeNormals
            (
                mesh,
                pp,
                patchEdges,
                coupledEdges
            )
        );

        forAll(en, patchEdgeI)
        {
            const edge& patchE = pp.edges()[patchEdgeI];
            // str.write(pp.localPoints()[pointi], en[pointi]);
            const point pt = patchE.centre(pp.localPoints());
            str.write(linePointRef(pt, pt + 0.1*en[patchEdgeI]));
        }
    }


    return 0;


//    {
//        OBJstream str(runTime.path()/"unweightedPatchNormals.obj");
//
//        const pointField pn
//        (
//            PatchTools::pointNormals
//            (
//                mesh,
//                pp,
//                identity(pp.size())+pp.start()
//            )
//        );
//        forAll(pn, pointi)
//        {
//            str.write(linePointRef(pp.localPoints()[pointi], pn[pointi]));
//        }
//    }
//    {
//        OBJstream str(runTime.path()/"areaWeightedPatchNormals.obj");
//
//        const pointField pn
//        (
//            areaPointNormals
//            (
//                mesh,
//                pp,
//                identity(pp.size())+pp.start()
//            )
//        );
//        forAll(pn, pointi)
//        {
//            str.write(linePointRef(pp.localPoints()[pointi], pn[pointi]));
//        }
//    }


    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
