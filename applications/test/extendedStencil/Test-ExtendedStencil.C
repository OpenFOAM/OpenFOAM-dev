/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    testExtendedStencil

Description
    Test app for determining extended stencil.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"
//#include "mapDistribute.H"
#include "OFstream.H"
#include "meshTools.H"
//#include "FECCellToFaceStencil.H"
//#include "CFCCellToFaceStencil.H"
//#include "CPCCellToFaceStencil.H"
//#include "CECCellToFaceStencil.H"
//#include "extendedCentredCellToFaceStencil.H"
//#include "extendedUpwindCellToFaceStencil.H"

//#include "centredCFCCellToFaceStencilObject.H"
//#include "centredFECCellToFaceStencilObject.H"
//#include "centredCPCCellToFaceStencilObject.H"
//#include "centredCECCellToFaceStencilObject.H"

//#include "upwindFECCellToFaceStencilObject.H"
//#include "upwindCPCCellToFaceStencilObject.H"
//#include "upwindCECCellToFaceStencilObject.H"

//#include "upwindCFCCellToFaceStencilObject.H"
//#include "centredCFCFaceToCellStencilObject.H"

#include "centredCECCellToCellStencilObject.H"
#include "centredCFCCellToCellStencilObject.H"
#include "centredCPCCellToCellStencilObject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeStencilOBJ
(
    const fileName& fName,
    const point& fc,
    const List<point>& stencilCc
)
{
    OFstream str(fName);
    label vertI = 0;

    meshTools::writeOBJ(str, fc);
    vertI++;

    forAll(stencilCc, i)
    {
        meshTools::writeOBJ(str, stencilCc[i]);
        vertI++;
        str << "l 1 " << vertI << nl;
    }
}


// Stats
void writeStencilStats(const labelListList& stencil)
{
    label sumSize = 0;
    label nSum = 0;
    label minSize = labelMax;
    label maxSize = labelMin;

    forAll(stencil, i)
    {
        const labelList& sCells = stencil[i];

        if (sCells.size() > 0)
        {
            sumSize += sCells.size();
            nSum++;
            minSize = min(minSize, sCells.size());
            maxSize = max(maxSize, sCells.size());
        }
    }
    reduce(sumSize, sumOp<label>());
    reduce(nSum, sumOp<label>());
    sumSize /= nSum;

    reduce(minSize, minOp<label>());
    reduce(maxSize, maxOp<label>());

    Info<< "Stencil size :" << nl
        << "    average : " << sumSize << nl
        << "    min     : " << minSize << nl
        << "    max     : " << maxSize << nl
        << endl;
}


// Main program:

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // Get times list
    instantList Times = runTime.times();
    #include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);
    #include "createMesh.H"

    // Force calculation of extended edge addressing
    const labelListList& edgeFaces = mesh.edgeFaces();
    const labelListList& edgeCells = mesh.edgeCells();
    const labelListList& pointCells = mesh.pointCells();
    Info<< "dummy:" << edgeFaces.size() + edgeCells.size() + pointCells.size()
        << endl;


    // Centred, semi-extended stencil (edge cells only)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//    {
//        // const FECCellToFaceStencil cfcStencil(mesh);
//        // const extendedCentredCellToFaceStencil addressing
//        //(
//        //    cfcStencil
//        //);
//        const extendedCentredStencil& addressing =
//        centredFECCellToFaceStencilObject::New
//        (
//            mesh
//        );
//
//        Info<< "faceEdgeCell:" << endl;
//        writeStencilStats(addressing.stencil());
//
//        // Collect stencil cell centres
//        List<List<point>> stencilPoints(mesh.nFaces());
//        addressing.collectData
//        (
//            mesh.C(),
//            stencilPoints
//        );
//
//        forAll(stencilPoints, facei)
//        {
//            writeStencilOBJ
//            (
//                runTime.path()/"faceEdgeCell" + Foam::name(facei) + ".obj",
//                mesh.faceCentres()[facei],
//                stencilPoints[facei]
//            );
//        }
//    }




//    // Centred, face stencil
//    // ~~~~~~~~~~~~~~~~~~~~~
//
//    {
//        const extendedCentredCellToFaceStencil& addressing =
//        centredCFCCellToFaceStencilObject::New
//        (
//            mesh
//        );
//
//        Info<< "cellFaceCell:" << endl;
//        writeStencilStats(addressing.stencil());
//
//
//        //// Do some interpolation.
//        //{
//        //    const labelListList& stencil = addressing.stencil();
//        //    List<List<scalar>> stencilWeights(stencil.size());
//        //    forAll(stencil, facei)
//        //    {
//        //        const labelList& fStencil = stencil[facei];
//        //
//        //        if (fStencil.size() > 0)
//        //        {
//        //            // Uniform weights
//        //            stencilWeights[facei] = scalarList
//        //            (
//        //                fStencil.size(),
//        //                1.0/fStencil.size()
//        //            );
//        //        }
//        //    }
//        //
//        //    tmp<surfaceVectorField> tfc
//        //    (
//        //        addressing.weightedSum(mesh.C(), stencilWeights)
//        //    );
//        //}
//
//
//        // Collect stencil cell centres
//        List<List<point>> stencilPoints(mesh.nFaces());
//        addressing.collectData
//        (
//            mesh.C(),
//            stencilPoints
//        );
//
//        forAll(stencilPoints, facei)
//        {
//            if (stencilPoints[facei].size() >= 15)
//            {
//                writeStencilOBJ
//                (
//                    runTime.path()/"centredFace" + Foam::name(facei) + ".obj",
//                    mesh.faceCentres()[facei],
//                    stencilPoints[facei]
//                );
//            }
//        }
//    }


//    // Centred, point stencil
//    // ~~~~~~~~~~~~~~~~~~~~~~
//
//    {
//        // const extendedCentredCellToFaceStencil& addressing =
//        // centredCPCStencilObject::New
//        //(
//        //    mesh
//        //);
//        //
//        // Info<< "cellPointCell:" << endl;
//        // writeStencilStats(addressing.stencil());
//        //
//        //
//        //// Collect stencil cell centres
//        // List<List<point>> stencilPoints(mesh.nFaces());
//        // addressing.collectData
//        //(
//        //    mesh.C(),
//        //    stencilPoints
//        //);
//        //
//        // forAll(stencilPoints, facei)
//        //{
//        //    writeStencilOBJ
//        //    (
//        //        runTime.path()/"centredPoint" + Foam::name(facei) + ".obj",
//        //        mesh.faceCentres()[facei],
//        //        stencilPoints[facei]
//        //    );
//        //}
//    }



//    // Centred, edge stencil
//    // ~~~~~~~~~~~~~~~~~~~~~~
//
//    {
//        // const extendedCentredCellToFaceStencil& addressing =
//        // centredCECStencilObject::New
//        //(
//        //    mesh
//        //);
//        //
//        // Info<< "cellEdgeCell:" << endl;
//        // writeStencilStats(addressing.stencil());
//        //
//        //
//        //// Collect stencil cell centres
//        // List<List<point>> stencilPoints(mesh.nFaces());
//        // addressing.collectData
//        //(
//        //    mesh.C(),
//        //    stencilPoints
//        //);
//        //
//        // forAll(stencilPoints, facei)
//        //{
//        //    writeStencilOBJ
//        //    (
//        //        runTime.path()/"centredEdge" + Foam::name(facei) + ".obj",
//        //        mesh.faceCentres()[facei],
//        //        stencilPoints[facei]
//        //    );
//        //}
//    }



    // Upwind, semi-extended stencil
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //{
    //    const extendedUpwindCellToFaceStencil& addressing =
    //    upwindFECCellToFaceStencilObject::New
    //    (
    //        mesh,
    //        0.5
    //    );
    //
    //    Info<< "upwind-faceEdgeCell:" << endl;
    //    writeStencilStats(addressing.ownStencil());
    //
    //    {
    //        // Collect stencil cell centres
    //        List<List<point>> ownPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.ownMap(),
    //            addressing.ownStencil(),
    //            mesh.C(),
    //            ownPoints
    //        );
    //
    //        forAll(ownPoints, facei)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"ownFEC" + Foam::name(facei) + ".obj",
    //                mesh.faceCentres()[facei],
    //                ownPoints[facei]
    //            );
    //        }
    //    }
    //    {
    //        // Collect stencil cell centres
    //        List<List<point>> neiPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.neiMap(),
    //            addressing.neiStencil(),
    //            mesh.C(),
    //            neiPoints
    //        );
    //
    //        forAll(neiPoints, facei)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"neiFEC" + Foam::name(facei) + ".obj",
    //                mesh.faceCentres()[facei],
    //                neiPoints[facei]
    //            );
    //        }
    //    }
    //}



    // Upwind, extended stencil
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    //{
    //    const extendedUpwindCellToFaceStencil& addressing =
    //    upwindCFCCellToFaceStencilObject::New
    //    (
    //        mesh,
    //        0.5
    //    );
    //
    //    Info<< "upwind-cellFaceCell:" << endl;
    //    writeStencilStats(addressing.ownStencil());
    //
    //    {
    //        // Collect stencil cell centres
    //        List<List<point>> ownPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.ownMap(),
    //            addressing.ownStencil(),
    //            mesh.C(),
    //            ownPoints
    //        );
    //
    //        forAll(ownPoints, facei)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"ownCFC" + Foam::name(facei) + ".obj",
    //                mesh.faceCentres()[facei],
    //                ownPoints[facei]
    //            );
    //        }
    //    }
    //    {
    //        // Collect stencil cell centres
    //        List<List<point>> neiPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.neiMap(),
    //            addressing.neiStencil(),
    //            mesh.C(),
    //            neiPoints
    //        );
    //
    //        forAll(neiPoints, facei)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"neiCFC" + Foam::name(facei) + ".obj",
    //                mesh.faceCentres()[facei],
    //                neiPoints[facei]
    //            );
    //        }
    //    }
    //}



    //---- CELL CENTRED STENCIL -----

    // Centred, cell stencil
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        const extendedCentredCellToCellStencil& addressing =
        centredCECCellToCellStencilObject::New
        (
            mesh
        );

        Info<< "cellCellCell:" << endl;
        writeStencilStats(addressing.stencil());

        // Collect stencil cell centres
        List<List<point>> stencilPoints(mesh.nCells());
        addressing.collectData
        (
            mesh.C(),
            stencilPoints
        );

        forAll(stencilPoints, celli)
        {
            writeStencilOBJ
            (
                runTime.path()/"centredCECCell" + Foam::name(celli) + ".obj",
                mesh.cellCentres()[celli],
                stencilPoints[celli]
            );
        }
    }
    {
        const extendedCentredCellToCellStencil& addressing =
        centredCPCCellToCellStencilObject::New
        (
            mesh
        );

        Info<< "cellCellCell:" << endl;
        writeStencilStats(addressing.stencil());

        // Collect stencil cell centres
        List<List<point>> stencilPoints(mesh.nCells());
        addressing.collectData
        (
            mesh.C(),
            stencilPoints
        );

        forAll(stencilPoints, celli)
        {
            writeStencilOBJ
            (
                runTime.path()/"centredCPCCell" + Foam::name(celli) + ".obj",
                mesh.cellCentres()[celli],
                stencilPoints[celli]
            );
        }
    }
    {
        const extendedCentredCellToCellStencil& addressing =
        centredCFCCellToCellStencilObject::New
        (
            mesh
        );

        Info<< "cellCellCell:" << endl;
        writeStencilStats(addressing.stencil());

        // Collect stencil cell centres
        List<List<point>> stencilPoints(mesh.nCells());
        addressing.collectData
        (
            mesh.C(),
            stencilPoints
        );

        forAll(stencilPoints, celli)
        {
            writeStencilOBJ
            (
                runTime.path()/"centredCFCCell" + Foam::name(celli) + ".obj",
                mesh.cellCentres()[celli],
                stencilPoints[celli]
            );
        }
    }


//XXXXXX
//    // Evaluate
//    List<List<scalar>> stencilData(faceStencils.size());
//    collectStencilData
//    (
//        distMap,
//        faceStencils,
//        vf,
//        stencilData
//    );
//    for (label faci = 0; faci < mesh.nInternalFaces(); faci++)
//    {
//        const scalarList& stData = stencilData[facei];
//        const scalarList& stWeight = fit[facei];
//
//        forAll(stData, i)
//        {
//            sf[facei] += stWeight[i]*stData[i];
//        }
//    }
//    See finiteVolume/lnInclude/leastSquaresGrad.C
//XXXXXX

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
