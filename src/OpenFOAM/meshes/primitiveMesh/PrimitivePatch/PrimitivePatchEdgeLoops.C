/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Create the list of loops of outside vertices. Goes wrong on multiply
    connected edges (loops will be unclosed).

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcEdgeLoops() const
{
    if (debug)
    {
        InfoInFunction << "Calculating boundary edge loops" << endl;
    }

    if (edgeLoopsPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "edge loops already calculated"
            << abort(FatalError);
    }

    const edgeList& patchEdges = edges();
    label nIntEdges = nInternalEdges();
    label nBdryEdges = patchEdges.size() - nIntEdges;

    if (nBdryEdges == 0)
    {
        edgeLoopsPtr_ = new labelListList(0);
        return;
    }

    const labelListList& patchPointEdges = pointEdges();


    //
    // Walk point-edge-point and assign loop number
    //

    // Loop per (boundary) edge.
    labelList loopNumber(nBdryEdges, -1);

    // Size return list plenty big
    edgeLoopsPtr_ = new labelListList(nBdryEdges);
    labelListList& edgeLoops = *edgeLoopsPtr_;


    // Current loop number.
    label loopI = 0;

    while (true)
    {
        // Find edge not yet given a loop number.
        label currentEdgeI = -1;

        for (label edgeI = nIntEdges; edgeI < patchEdges.size(); edgeI++)
        {
            if (loopNumber[edgeI-nIntEdges] == -1)
            {
                currentEdgeI = edgeI;
                break;
            }
        }

        if (currentEdgeI == -1)
        {
            // Did not find edge not yet assigned a loop number so done all.
            break;
        }

        // Temporary storage for vertices of current loop
        DynamicList<label> loop(nBdryEdges);

        // Walk from first all the way round, assigning loops
        label currentVertI = patchEdges[currentEdgeI].start();

        do
        {
            loop.append(currentVertI);

            loopNumber[currentEdgeI - nIntEdges] = loopI;

            // Step to next vertex
            currentVertI = patchEdges[currentEdgeI].otherVertex(currentVertI);

            // Step to next (unmarked, boundary) edge.
            const labelList& curEdges = patchPointEdges[currentVertI];

            currentEdgeI = -1;

            forAll(curEdges, pI)
            {
                label edgeI = curEdges[pI];

                if (edgeI >= nIntEdges && (loopNumber[edgeI - nIntEdges] == -1))
                {
                    // Unassigned boundary edge.
                    currentEdgeI = edgeI;

                    break;
                }
            }
        }
        while (currentEdgeI != -1);

        // Done all for current loop. Transfer to edgeLoops.
        edgeLoops[loopI].transfer(loop);

        loopI++;
    }

    edgeLoops.setSize(loopI);

    if (debug)
    {
        Info<< "    Finished." << endl;
    }
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::edgeLoops() const
{
    if (!edgeLoopsPtr_)
    {
        calcEdgeLoops();
    }

    return *edgeLoopsPtr_;
}


// ************************************************************************* //
