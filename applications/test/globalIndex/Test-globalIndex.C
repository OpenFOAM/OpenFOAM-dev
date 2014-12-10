/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    globalIndexTest

Description
    Simple demonstration and test application for the globalIndex class.

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IOstreams.H"
#include "OStringStream.H"
#include "IStringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    // Global numbering of cells (proc0 elements first, then proc1, etc.)
    globalIndex globalNumbering(mesh.nCells());

    if (globalNumbering.localSize() != mesh.nCells())
    {
        FatalErrorIn(args.executable())
            << "Problem." << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        WarningIn(args.executable())
            << "globalIndex class is only useful in parallel code."
            << endl;
    }

    // convert from local to global and back.
    for (label cellI = 0; cellI < mesh.nCells(); cellI++)
    {
        // to global index
        label globalCellI = globalNumbering.toGlobal(cellI);

        // and back
        label procI = globalNumbering.whichProcID(globalCellI);
        label localCellI = globalNumbering.toLocal(globalCellI);

        if (procI != Pstream::myProcNo() || localCellI != cellI)
        {
            FatalErrorIn(args.executable())
                << "Problem. cellI:" << cellI << " localCellI:" << localCellI
                << " procI:" << procI << abort(FatalError);
        }

        if (!globalNumbering.isLocal(globalCellI))
        {
            FatalErrorIn(args.executable())
                << "Problem. cellI:" << cellI << " globalCellI:" << globalCellI
                << " not local" << abort(FatalError);
        }
    }


    // Try whichProcID on a few borderline cases.

    if (mesh.nCells() < 1)
    {
        FatalErrorIn(args.executable())
            << "Test needs to be run on a case with at least one"
            << " cell per processor." << abort(FatalError);
    }

    if (Pstream::myProcNo() > 0)
    {
        // We already checked that toGlobal(0) maps back correctly to myProcNo
        // so now check that the index one before maps to the previous processor
        label prevProcCellI = globalNumbering.toGlobal(0)-1;
        label procI = globalNumbering.whichProcID(prevProcCellI);

        if (procI != Pstream::myProcNo()-1)
        {
            FatalErrorIn(args.executable())
                << "Problem. global:" << prevProcCellI
                << " expected on processor:" << Pstream::myProcNo()-1
                << " but is calculated to be on procI:" << procI
                << abort(FatalError);
        }

        if (globalNumbering.isLocal(prevProcCellI))
        {
            FatalErrorIn(args.executable())
                << "Problem. globalCellI:" << prevProcCellI
                << " calculated as local" << abort(FatalError);
        }

        if (!globalNumbering.isLocal(procI, prevProcCellI))
        {
            FatalErrorIn(args.executable())
                << "Problem. globalCellI:" << prevProcCellI
                << " not calculated as local on processor:" << procI
                << abort(FatalError);
        }
    }


    if (Pstream::myProcNo() < Pstream::nProcs()-1)
    {
        label nextProcCellI = globalNumbering.toGlobal(mesh.nCells()-1)+1;
        label procI = globalNumbering.whichProcID(nextProcCellI);

        if (procI != Pstream::myProcNo()+1)
        {
            FatalErrorIn(args.executable())
                << "Problem. global:" << nextProcCellI
                << " expected on processor:" << Pstream::myProcNo()+1
                << " but is calculated to be on procI:" << procI
                << abort(FatalError);
        }

        if (globalNumbering.isLocal(nextProcCellI))
        {
            FatalErrorIn(args.executable())
                << "Problem. globalCellI:" << nextProcCellI
                << " calculated as local" << abort(FatalError);
        }

        if (!globalNumbering.isLocal(procI, nextProcCellI))
        {
            FatalErrorIn(args.executable())
                << "Problem. globalCellI:" << nextProcCellI
                << " not calculated as local on processor:" << procI
                << abort(FatalError);
        }
    }

    return 0;
}


// ************************************************************************* //
