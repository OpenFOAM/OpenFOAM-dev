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

\*---------------------------------------------------------------------------*/

#include "processorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, processorFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::processorFvPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        const vectorField delta(coupledFvPatch::delta());

        // The face normals point in the opposite direction on the other side
        const vectorField nbrDelta
        (
            procPolyPatch_.neighbFaceCentres()
          - procPolyPatch_.neighbFaceCellCentres()
        );

        const scalarField nfDelta(nf() & delta);

        const scalarField nbrNfDelta
        (
            (
                procPolyPatch_.neighbFaceAreas()
               /(mag(procPolyPatch_.neighbFaceAreas()) + vSmall)
            ) & nbrDelta
        );

        forAll(delta, facei)
        {
            const scalar ndoi = nfDelta[facei];
            const scalar ndni = nbrNfDelta[facei];
            const scalar ndi = ndoi + ndni;

            if (ndni/vGreat < ndi)
            {
                w[facei] = ndni/ndi;
            }
            else
            {
                const scalar doi = mag(delta[facei]);
                const scalar dni = mag(nbrDelta[facei]);
                const scalar di = doi + dni;

                w[facei] = dni/di;
            }
        }
    }
    else
    {
        w = 1.0;
    }
}


Foam::tmp<Foam::vectorField> Foam::processorFvPatch::delta() const
{
    if (Pstream::parRun())
    {
        // To the transformation if necessary
        if (parallel())
        {
            return
                coupledFvPatch::delta()
              - (
                    procPolyPatch_.neighbFaceCentres()
                  - procPolyPatch_.neighbFaceCellCentres()
                );
        }
        else
        {
            return
                coupledFvPatch::delta()
              - transform
                (
                    forwardT(),
                    (
                        procPolyPatch_.neighbFaceCentres()
                      - procPolyPatch_.neighbFaceCellCentres()
                    )
                );
        }
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::processorFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::processorFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}


// ************************************************************************* //
