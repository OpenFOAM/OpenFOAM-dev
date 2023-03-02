/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "cutPolyValue.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
Type min(const UIndirectList<Type>& l)
{
    Type result = pTraits<Type>::max;
    forAll(l, i)
    {
        result = min(result, l[i]);
    }
    return result;
}

template<class Type>
Type max(const UIndirectList<Type>& l)
{
    Type result = pTraits<Type>::min;
    forAll(l, i)
    {
        result = max(result, l[i]);
    }
    return result;
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::labelPair> Foam::cutPoly::faceCuts
(
    const face& f,
    const scalarField& pAlphas,
    const scalar isoAlpha
)
{
    UIndirectList<scalar> fAlphas(pAlphas, f);

    // Quick reject if all alpha values are above or below the iso value
    if (min(fAlphas) > isoAlpha || isoAlpha > max(fAlphas))
    {
        return List<labelPair>();
    }

    // Find the starting point
    label fpi0 = 0;
    while ((fAlphas[fpi0] < isoAlpha) == separatedBelow)
    {
        ++ fpi0;
    }

    // Create cuts on every edge for which there is a sign change
    DynamicList<labelPair> cuts(1);
    label cuti = 0;
    forAll(f, i)
    {
        const label fpi1 = f.fcIndex(fpi0);

        if ((fAlphas[fpi0] < isoAlpha) != (fAlphas[fpi1] < isoAlpha))
        {
            if (cuti == 0)
            {
                cuts.append({fpi0, -1});
            }
            else
            {
                cuts.last()[1] = fpi0;
            }

            cuti = cuts.last().fcIndex(cuti);
        }

        fpi0 = fpi1;
    }

    // There should have been an even number of cuts
    if (cuti != 0)
    {
        FatalErrorInFunction
            << "Cutting values " << fAlphas << " with iso-value " << isoAlpha
            << " resulted in an odd number of cuts"
            << exit(FatalError);
    }

    cuts.shrink();
    return cuts;
}


Foam::labelListList Foam::cutPoly::cellCuts
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const faceUList& fs,
    const List<List<labelPair>>& fCuts,
    const scalarField& pAlphas,
    const scalar isoAlpha
)
{
    // Quick return if not cut
    label nCellFaceCuts = 0;
    forAll(c, cfi)
    {
        nCellFaceCuts += fCuts[c[cfi]].size();
    }
    if (nCellFaceCuts == 0)
    {
        return labelListList();
    }

    // Get local cell-face-edge addressing
    const CompactListList<label>& cfiAndFeiToCei = cAddr.cfiAndFeiToCei();
    const List<Pair<labelPair>>& ceiToCfiAndFei = cAddr.ceiToCfiAndFei();
    const boolList& cOwns = cAddr.cOwns();

    // For each cell-face and face-edge, get the face-edge at the other end of
    // the edge's cut, or -1 if the edge is not cut. This allows us to walk
    // along face cuts. It also doubles as the "visited" array during the walk.
    // As cuts are visited, the relevant label in this array is set to -1 to
    // ensure that the cut is not considered twice.
    labelListList cellFaceAndFaceEdgeCutToFaceEdge(c.size());
    forAll(c, cfi)
    {
        cellFaceAndFaceEdgeCutToFaceEdge[cfi] =
            labelList(fs[c[cfi]].size(), -1);

        forAll(fCuts[c[cfi]], fci)
        {
            const label fei0 = fCuts[c[cfi]][fci].first();
            const label fei1 = fCuts[c[cfi]][fci].second();

            if (cOwns[cfi])
            {
                cellFaceAndFaceEdgeCutToFaceEdge[cfi][fei0] = fei1;
            }
            else
            {
                cellFaceAndFaceEdgeCutToFaceEdge[cfi][fei1] = fei0;
            }
        }
    }

    // Walk around the face cuts to generate the cell cuts
    DynamicList<labelList> cuts;
    {
        label nCellFaceCutsVisited = 0;
        label cfi0 = 0, fei0 = 0;
        while (nCellFaceCutsVisited < nCellFaceCuts)
        {
            // Find the next unvisited connection
            bool found = false;
            for (label cfj0 = cfi0; cfj0 < c.size(); ++ cfj0)
            {
                for (label fej0 = fei0; fej0 < fs[c[cfj0]].size(); ++ fej0)
                {
                    const label fej1 =
                        cellFaceAndFaceEdgeCutToFaceEdge[cfj0][fej0];

                    if (fej1 != -1)
                    {
                        cfi0 = cfj0;
                        fei0 = fej0;
                        found = true;
                    }

                    if (found) break;
                }

                if (found) break;

                fei0 = 0;
            }
            if (!found)
            {
                FatalErrorInFunction
                    << "Could not find next unvisited connection for cell cut"
                    << exit(FatalError);
            }

            // Walk around the cell from face to face to form a cut
            label cfi = cfi0, fei = fei0;
            DynamicList<label> cut(8);
            while (cellFaceAndFaceEdgeCutToFaceEdge[cfi][fei] != -1)
            {
                ++ nCellFaceCutsVisited;

                const label otherFei =
                    cellFaceAndFaceEdgeCutToFaceEdge[cfi][fei];

                cellFaceAndFaceEdgeCutToFaceEdge[cfi][fei] = -1;

                const label cei = cfiAndFeiToCei[cfi][otherFei];

                cut.append(cei);

                const labelPair next =
                    ceiToCfiAndFei[cei].other({cfi, otherFei});

                cfi = next.first();
                fei = next.second();
            }

            // Copy the data over into the result list
            cuts.append(labelList());
            cuts.last().transfer(cut);
        }
    }

    cuts.shrink();
    return cuts;
}


void Foam::cutPoly::writeFaceCuts
(
    const face& f,
    const List<labelPair>& fCuts,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    OBJstream& obj
)
{
    forAll(fCuts, cuti)
    {
        pointField cutPs(2);

        forAll(fCuts[cuti], i)
        {
            const edge e = f.faceEdge(fCuts[cuti][i]);

            cutPs[i] = edgeCutValue(e, pAlphas, isoAlpha, ps);
        }

        obj.write(edge(0, 1), cutPs);
    }
}


void Foam::cutPoly::writeCellCuts
(
    const cell& c,
    const cellEdgeAddressing& cAddr,
    const List<List<label>>& cCuts,
    const faceList& fs,
    const pointField& ps,
    const scalarField& pAlphas,
    const scalar isoAlpha,
    OBJstream& obj
)
{
    // Quick return if not cut
    if (cCuts.size() == 0)
    {
        return;
    }

    // Get local cell-face-edge addressing
    const List<Pair<labelPair>>& ceiToCfiAndFei = cAddr.ceiToCfiAndFei();

    // Write out each cut as a single face
    forAll(cCuts, cuti)
    {
        if (cCuts[cuti].size() < 3) continue;

        pointField cutPs(cCuts[cuti].size());

        forAll(cCuts[cuti], i)
        {
            const label cei = cCuts[cuti][i];
            const label cfi = ceiToCfiAndFei[cei][0][0];
            const label fei = ceiToCfiAndFei[cei][0][1];

            const edge e = fs[c[cfi]].faceEdge(fei);

            cutPs[i] = edgeCutValue(e, pAlphas, isoAlpha, ps);
        }

        obj.write(face(identityMap(cCuts[cuti].size())), cutPs, false);
    }
}


// ************************************************************************* //
