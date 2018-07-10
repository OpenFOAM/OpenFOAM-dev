/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    // Global numbering of cells (proc0 elements first, then proc1, etc.)
    globalIndex globalNumbering(mesh.nCells());

    if (globalNumbering.localSize() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Problem." << abort(FatalError);
    }


    if (!Pstream::parRun())
    {
        WarningInFunction
            << "globalIndex class is only useful in parallel code."
            << endl;
    }

    // convert from local to global and back.
    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        // to global index
        label globalCelli = globalNumbering.toGlobal(celli);

        // and back
        label proci = globalNumbering.whichProcID(globalCelli);
        label localCelli = globalNumbering.toLocal(globalCelli);

        if (proci != Pstream::myProcNo() || localCelli != celli)
        {
            FatalErrorInFunction
                << "Problem. celli:" << celli << " localCelli:" << localCelli
                << " proci:" << proci << abort(FatalError);
        }

        if (!globalNumbering.isLocal(globalCelli))
        {
            FatalErrorInFunction
                << "Problem. celli:" << celli << " globalCelli:" << globalCelli
                << " not local" << abort(FatalError);
        }
    }


    // Try whichProcID on a few borderline cases.

    if (mesh.nCells() < 1)
    {
        FatalErrorInFunction
            << "Test needs to be run on a case with at least one"
            << " cell per processor." << abort(FatalError);
    }

    if (Pstream::myProcNo() > 0)
    {
        // We already checked that toGlobal(0) maps back correctly to myProcNo
        // so now check that the index one before maps to the previous processor
        label prevProcCelli = globalNumbering.toGlobal(0)-1;
        label proci = globalNumbering.whichProcID(prevProcCelli);

        if (proci != Pstream::myProcNo()-1)
        {
            FatalErrorInFunction
                << "Problem. global:" << prevProcCelli
                << " expected on processor:" << Pstream::myProcNo()-1
                << " but is calculated to be on proci:" << proci
                << abort(FatalError);
        }

        if (globalNumbering.isLocal(prevProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << prevProcCelli
                << " calculated as local" << abort(FatalError);
        }

        if (!globalNumbering.isLocal(proci, prevProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << prevProcCelli
                << " not calculated as local on processor:" << proci
                << abort(FatalError);
        }
    }


    if (Pstream::myProcNo() < Pstream::nProcs()-1)
    {
        label nextProcCelli = globalNumbering.toGlobal(mesh.nCells()-1)+1;
        label proci = globalNumbering.whichProcID(nextProcCelli);

        if (proci != Pstream::myProcNo()+1)
        {
            FatalErrorInFunction
                << "Problem. global:" << nextProcCelli
                << " expected on processor:" << Pstream::myProcNo()+1
                << " but is calculated to be on proci:" << proci
                << abort(FatalError);
        }

        if (globalNumbering.isLocal(nextProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << nextProcCelli
                << " calculated as local" << abort(FatalError);
        }

        if (!globalNumbering.isLocal(proci, nextProcCelli))
        {
            FatalErrorInFunction
                << "Problem. globalCelli:" << nextProcCelli
                << " not calculated as local on processor:" << proci
                << abort(FatalError);
        }
    }

    return 0;
}


// ************************************************************************* //
