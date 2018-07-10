/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "cyclicACMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicACMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cyclicACMIFvPatch::updateAreas() const
{
    if (cyclicACMIPatch().updated())
    {
        if (debug)
        {
            Pout<< "cyclicACMIFvPatch::updateAreas() : updating fv areas for "
                << name() << " and " << nonOverlapFvPatch().name()
                << endl;
        }

        // owner couple
        const_cast<vectorField&>(Sf()) = patch().faceAreas();
        const_cast<scalarField&>(magSf()) = mag(patch().faceAreas());

        // owner non-overlapping
        const fvPatch& nonOverlapPatch = nonOverlapFvPatch();
        const_cast<vectorField&>(nonOverlapPatch.Sf()) =
            nonOverlapPatch.patch().faceAreas();
        const_cast<scalarField&>(nonOverlapPatch.magSf()) =
            mag(nonOverlapPatch.patch().faceAreas());

        // neighbour couple
        const cyclicACMIFvPatch& nbrACMI = neighbFvPatch();
        const_cast<vectorField&>(nbrACMI.Sf()) =
            nbrACMI.patch().faceAreas();
        const_cast<scalarField&>(nbrACMI.magSf()) =
            mag(nbrACMI.patch().faceAreas());

        // neighbour non-overlapping
        const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapFvPatch();
        const_cast<vectorField&>(nbrNonOverlapPatch.Sf()) =
            nbrNonOverlapPatch.patch().faceAreas();
        const_cast<scalarField&>(nbrNonOverlapPatch.magSf()) =
            mag(nbrNonOverlapPatch.patch().faceAreas());

        // set the updated flag
        cyclicACMIPatch().setUpdated(false);
    }
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicACMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        // These deltas are of the cyclic part alone - they are
        // not affected by the amount of overlap with the nonOverlapPatch
        const scalarField deltan(this->deltan());
        const scalarField nbrDeltan(this->nbrDeltan());

        forAll(deltan, facei)
        {
            scalar di = deltan[facei];
            scalar dni = nbrDeltan[facei];

            if (dni < cyclicACMIPolyPatch::tolerance())
            {
                // Avoid zero weights on disconnected faces. This value
                // will be weighted with the (zero) face area so will not
                // influence calculations.
                w[facei] = 1.0;
            }
            else
            {
                w[facei] = dni/(di + dni);
            }
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


// ************************************************************************* //
