/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "cyclicAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicAMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicAMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


void Foam::cyclicAMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const cyclicAMIFvPatch& nbrPatch = neighbFvPatch();

        const scalarField deltas(nf() & coupledFvPatch::delta());

        tmp<scalarField> tnbrDeltas;
        if (applyLowWeightCorrection())
        {
            tnbrDeltas =
                interpolate
                (
                    nbrPatch.nf() & nbrPatch.coupledFvPatch::delta(),
                    scalarField(this->size(), 1.0)
                );
        }
        else
        {
            tnbrDeltas =
                interpolate(nbrPatch.nf() & nbrPatch.coupledFvPatch::delta());
        }

        const scalarField& nbrDeltas = tnbrDeltas();

        forAll(deltas, faceI)
        {
            scalar di = deltas[faceI];
            scalar dni = nbrDeltas[faceI];

            w[faceI] = dni/(di + dni);
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


Foam::tmp<Foam::vectorField> Foam::cyclicAMIFvPatch::delta() const
{
    const cyclicAMIFvPatch& nbrPatch = neighbFvPatch();

    if (coupled())
    {
        const vectorField patchD(coupledFvPatch::delta());

        tmp<vectorField> tnbrPatchD;
        if (applyLowWeightCorrection())
        {
            tnbrPatchD =
                interpolate
                (
                    nbrPatch.coupledFvPatch::delta(),
                    vectorField(this->size(), vector::zero)
                );
        }
        else
        {
            tnbrPatchD = interpolate(nbrPatch.coupledFvPatch::delta());
        }

        const vectorField& nbrPatchD = tnbrPatchD();

        tmp<vectorField> tpdv(new vectorField(patchD.size()));
        vectorField& pdv = tpdv();

        // do the transformation if necessary
        if (parallel())
        {
            forAll(patchD, faceI)
            {
                const vector& ddi = patchD[faceI];
                const vector& dni = nbrPatchD[faceI];

                pdv[faceI] = ddi - dni;
            }
        }
        else
        {
            forAll(patchD, faceI)
            {
                const vector& ddi = patchD[faceI];
                const vector& dni = nbrPatchD[faceI];

                pdv[faceI] = ddi - transform(forwardT()[0], dni);
            }
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicAMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// ************************************************************************* //
