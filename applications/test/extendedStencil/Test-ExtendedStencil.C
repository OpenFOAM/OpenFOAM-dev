/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();
#   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);
#   include "createMesh.H"

    // Force calculation of extended edge addressing
    const labelListList& edgeFaces = mesh.edgeFaces();
    const labelListList& edgeCells = mesh.edgeCells();
    const labelListList& pointCells = mesh.pointCells();
    Info<< "dummy:" << edgeFaces.size() + edgeCells.size() + pointCells.size()
        << endl;


    // Centred, semi-extended stencil (edge cells only)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//    {
//        //const FECCellToFaceStencil cfcStencil(mesh);
//        //const extendedCentredCellToFaceStencil addressing
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
//        List<List<point> > stencilPoints(mesh.nFaces());
//        addressing.collectData
//        (
//            mesh.C(),
//            stencilPoints
//        );
//
//        forAll(stencilPoints, faceI)
//        {
//            writeStencilOBJ
//            (
//                runTime.path()/"faceEdgeCell" + Foam::name(faceI) + ".obj",
//                mesh.faceCentres()[faceI],
//                stencilPoints[faceI]
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
//        //    List<List<scalar> > stencilWeights(stencil.size());
//        //    forAll(stencil, faceI)
//        //    {
//        //        const labelList& fStencil = stencil[faceI];
//        //
//        //        if (fStencil.size() > 0)
//        //        {
//        //            // Uniform weights
//        //            stencilWeights[faceI] = scalarList
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
//        List<List<point> > stencilPoints(mesh.nFaces());
//        addressing.collectData
//        (
//            mesh.C(),
//            stencilPoints
//        );
//
//        forAll(stencilPoints, faceI)
//        {
//            if (stencilPoints[faceI].size() >= 15)
//            {
//                writeStencilOBJ
//                (
//                    runTime.path()/"centredFace" + Foam::name(faceI) + ".obj",
//                    mesh.faceCentres()[faceI],
//                    stencilPoints[faceI]
//                );
//            }
//        }
//    }


//    // Centred, point stencil
//    // ~~~~~~~~~~~~~~~~~~~~~~
//
//    {
//        //const extendedCentredCellToFaceStencil& addressing =
//        //centredCPCStencilObject::New
//        //(
//        //    mesh
//        //);
//        //
//        //Info<< "cellPointCell:" << endl;
//        //writeStencilStats(addressing.stencil());
//        //
//        //
//        //// Collect stencil cell centres
//        //List<List<point> > stencilPoints(mesh.nFaces());
//        //addressing.collectData
//        //(
//        //    mesh.C(),
//        //    stencilPoints
//        //);
//        //
//        //forAll(stencilPoints, faceI)
//        //{
//        //    writeStencilOBJ
//        //    (
//        //        runTime.path()/"centredPoint" + Foam::name(faceI) + ".obj",
//        //        mesh.faceCentres()[faceI],
//        //        stencilPoints[faceI]
//        //    );
//        //}
//    }



//    // Centred, edge stencil
//    // ~~~~~~~~~~~~~~~~~~~~~~
//
//    {
//        //const extendedCentredCellToFaceStencil& addressing =
//        //centredCECStencilObject::New
//        //(
//        //    mesh
//        //);
//        //
//        //Info<< "cellEdgeCell:" << endl;
//        //writeStencilStats(addressing.stencil());
//        //
//        //
//        //// Collect stencil cell centres
//        //List<List<point> > stencilPoints(mesh.nFaces());
//        //addressing.collectData
//        //(
//        //    mesh.C(),
//        //    stencilPoints
//        //);
//        //
//        //forAll(stencilPoints, faceI)
//        //{
//        //    writeStencilOBJ
//        //    (
//        //        runTime.path()/"centredEdge" + Foam::name(faceI) + ".obj",
//        //        mesh.faceCentres()[faceI],
//        //        stencilPoints[faceI]
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
    //        List<List<point> > ownPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.ownMap(),
    //            addressing.ownStencil(),
    //            mesh.C(),
    //            ownPoints
    //        );
    //
    //        forAll(ownPoints, faceI)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"ownFEC" + Foam::name(faceI) + ".obj",
    //                mesh.faceCentres()[faceI],
    //                ownPoints[faceI]
    //            );
    //        }
    //    }
    //    {
    //        // Collect stencil cell centres
    //        List<List<point> > neiPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.neiMap(),
    //            addressing.neiStencil(),
    //            mesh.C(),
    //            neiPoints
    //        );
    //
    //        forAll(neiPoints, faceI)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"neiFEC" + Foam::name(faceI) + ".obj",
    //                mesh.faceCentres()[faceI],
    //                neiPoints[faceI]
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
    //        List<List<point> > ownPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.ownMap(),
    //            addressing.ownStencil(),
    //            mesh.C(),
    //            ownPoints
    //        );
    //
    //        forAll(ownPoints, faceI)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"ownCFC" + Foam::name(faceI) + ".obj",
    //                mesh.faceCentres()[faceI],
    //                ownPoints[faceI]
    //            );
    //        }
    //    }
    //    {
    //        // Collect stencil cell centres
    //        List<List<point> > neiPoints(mesh.nFaces());
    //        addressing.collectData
    //        (
    //            addressing.neiMap(),
    //            addressing.neiStencil(),
    //            mesh.C(),
    //            neiPoints
    //        );
    //
    //        forAll(neiPoints, faceI)
    //        {
    //            writeStencilOBJ
    //            (
    //                runTime.path()/"neiCFC" + Foam::name(faceI) + ".obj",
    //                mesh.faceCentres()[faceI],
    //                neiPoints[faceI]
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
        List<List<point> > stencilPoints(mesh.nCells());
        addressing.collectData
        (
            mesh.C(),
            stencilPoints
        );

        forAll(stencilPoints, cellI)
        {
            writeStencilOBJ
            (
                runTime.path()/"centredCECCell" + Foam::name(cellI) + ".obj",
                mesh.cellCentres()[cellI],
                stencilPoints[cellI]
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
        List<List<point> > stencilPoints(mesh.nCells());
        addressing.collectData
        (
            mesh.C(),
            stencilPoints
        );

        forAll(stencilPoints, cellI)
        {
            writeStencilOBJ
            (
                runTime.path()/"centredCPCCell" + Foam::name(cellI) + ".obj",
                mesh.cellCentres()[cellI],
                stencilPoints[cellI]
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
        List<List<point> > stencilPoints(mesh.nCells());
        addressing.collectData
        (
            mesh.C(),
            stencilPoints
        );

        forAll(stencilPoints, cellI)
        {
            writeStencilOBJ
            (
                runTime.path()/"centredCFCCell" + Foam::name(cellI) + ".obj",
                mesh.cellCentres()[cellI],
                stencilPoints[cellI]
            );
        }
    }


//XXXXXX
//    // Evaluate
//    List<List<scalar> > stencilData(faceStencils.size());
//    collectStencilData
//    (
//        distMap,
//        faceStencils,
//        vf,
//        stencilData
//    );
//    for (label faci = 0; faci < mesh.nInternalFaces(); faci++)
//    {
//        const scalarList& stData = stencilData[faceI];
//        const scalarList& stWeight = fit[faceI];
//
//        forAll(stData, i)
//        {
//            sf[faceI] += stWeight[i]*stData[i];
//        }
//    }
//    See finiteVolume/lnInclude/leastSquaresGrad.C
//XXXXXX

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
