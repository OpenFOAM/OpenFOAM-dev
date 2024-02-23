/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "argList.H"
#include "labelVector.H"
#include "ListOps.H"
#include "polyMesh.H"
#include "primitivePatchIntersection.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    static const scalar defaultSnapTol = Foam::sqrt(rootSmall);

    argList::validArgs.append("source");
    argList::validArgs.append("target");
    argList::addOption
    (
        "snapTol",
        "snapTol",
        "snapping tolerance - default is " + name(defaultSnapTol)
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const polyPatch& srcPatch = mesh.boundaryMesh()[args[1]];
    const polyPatch& tgtPatch = mesh.boundaryMesh()[args[2]];
    const scalar snapTol =
        args.optionLookupOrDefault<scalar>("snapTol", defaultSnapTol);

    if (mesh.nGeometricD() != 2 && mesh.nGeometricD() != 3)
    {
        FatalErrorInFunction
            << "Primitive patch intersection can only be done on 2 or "
            << "3-dimensional meshes" << exit(FatalError);
    }

    primitivePatchIntersection intersection(srcPatch, tgtPatch, snapTol);

    // Check the patch-edge-points start and end in the right places
    auto check = [&](const bool isSrc)
    {
        const primitivePatch& patch = isSrc ? srcPatch : tgtPatch;

        const labelList& patchPointPoints =
            isSrc
          ? intersection.srcPointPoints()
          : intersection.tgtPointPoints();

        const List<DynamicList<label>>& patchEdgePoints =
            isSrc
          ? intersection.srcEdgePoints()
          : intersection.tgtEdgePoints();

        forAll(patch.edges(), edgei)
        {
            const edge& patchE = patch.edges()[edgei];

            const edge ictE
            (
                patchPointPoints[patchE.start()],
                patchPointPoints[patchE.end()]
            );

            const labelList& ictEPs = patchEdgePoints[edgei];

            if (ictE != edge(ictEPs.first(), ictEPs.last()))
            {
                FatalErrorInFunction
                    << "Patch-edge " << patchE << " has intersection end "
                    << "points " << ictE << ", but these do not match the "
                    << "first and last of all the intersection points "
                    << ictEPs << exit(FatalError);
            }
        }
    };
    check(true);
    check(false);

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
