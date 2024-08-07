/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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

#include "symmetryPlanePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "RemoteData.H"
#include "Tuple3.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symmetryPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, dictionary);
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::symmetryPlanePolyPatch::calcGeometry(PstreamBuffers&)
{
    if (n_ != vector::rootMax)
    {
        return;
    }

    if (!returnReduce(size(), sumOp<label>()))
    {
        return;
    }

    const pointField& cf(primitivePatch::faceCentres());
    const vectorField& nf(faceNormals());
    n_ = gAverage(nf);

    if (debug)
    {
        Info<< "Patch " << name() << " calculated average normal "
            << n_ << endl;
    }

    // Get the largest variation from the average normal
    RemoteData<Tuple3<scalar, vector, vector>> maxDeltaN;
    if (size())
    {
        const scalarField deltaNSqr(magSqr(n_ - nf));
        const label i = findMax(deltaNSqr);
        maxDeltaN = {Pstream::myProcNo(), i, {deltaNSqr[i], cf[i], nf[i]}};
    }
    combineReduce
    (
        maxDeltaN,
        RemoteData<Tuple3<scalar, vector, vector>>::greatestFirstEqOp()
    );

    // Fail if the symmetry plane is not planar
    if (maxDeltaN.data.first() > small)
    {
        Ostream& fos =
            FatalErrorInFunction
                << "Symmetry plane '" << name()
                << "' is not planar" << endl;

        fos << "At patch face #" << maxDeltaN.elementi;

        if (Pstream::parRun())
        {
            fos << " on processor #" << maxDeltaN.proci;
        }

        fos << " with centre " << maxDeltaN.data.second()
            << " the normal " << maxDeltaN.data.third()
            << " differs from the average normal " << n_
            << " by " << maxDeltaN.data.first() << nl
            << "Either split the patch into planar parts or use the "
            << symmetryPolyPatch::typeName << " patch type"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    n_(vector::rootMax)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    n_(vector::rootMax)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    n_(pp.n_)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    n_(pp.n_)
{}


// ************************************************************************* //
