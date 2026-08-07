/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
#include "volFields.H"
#include "surfaceFields.H"

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
        if (!mesh().conformal())
        {
            coupledFvPatch::makeWeights
            (
                w,
              - mesh().Sf().boundaryField()[index()],
                mesh().Cf().boundaryField()[index()]
              - mesh().C().boundaryField()[index()]
            );
        }
        else
        {
            coupledFvPatch::makeWeights
            (
                w,
                processorPoly_.neighbFaceAreas(),
                processorPoly_.neighbFaceCentres()
              - processorPoly_.neighbFaceCellCentres()
            );
        }
    }
    else
    {
        w = 1;
    }
}


Foam::tmp<Foam::vectorField> Foam::processorFvPatch::delta() const
{
    if (Pstream::parRun())
    {
        if (!mesh().conformal())
        {
            return
                coupledFvPatch::delta()
              - (
                    mesh().Cf().boundaryField()[index()]
                  - mesh().C().boundaryField()[index()]
                );
        }
        else
        {
            return
                coupledFvPatch::delta()
              - transform().transform
                (
                    eval
                    (
                        processorPoly_.neighbFaceCentres()
                      - processorPoly_.neighbFaceCellCentres()
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
