/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "findRefCell.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::setRefCell
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    if (fieldRef.needReference() || forceReference)
    {
        word refCellName = field.name() + "RefCell";
        word refPointName = field.name() + "RefPoint";

        word refValueName = field.name() + "RefValue";

        if (dict.found(refCellName))
        {
            if (Pstream::master())
            {
                refCelli = dict.lookup<label>(refCellName);

                if (refCelli < 0 || refCelli >= field.mesh().nCells())
                {
                    FatalIOErrorInFunction
                    (
                        dict
                    )   << "Illegal master cellID " << refCelli
                        << ". Should be 0.." << field.mesh().nCells()
                        << exit(FatalIOError);
                }
            }
            else
            {
                refCelli = -1;
            }
        }
        else if (dict.found(refPointName))
        {
            point refPointi(dict.lookup(refPointName));

            // Try a linear search with approximate face-planes test to avoid
            // octree and tet-base-point construction
            refCelli =
                meshSearch::findCellNoTree
                (
                    field.mesh(),
                    refPointi,
                    pointInCellShapes::facePlanes
                );

            label hasRef = (refCelli >= 0 ? 1 : 0);
            label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

            // If a reference cell was not found then use a robust cell-tet
            // test with an octree search
            if (sumHasRef != 1)
            {
                refCelli = meshSearch::New(field.mesh()).findCell(refPointi);

                hasRef = (refCelli >= 0 ? 1 : 0);
                sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
            }

            if (sumHasRef != 1)
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Unable to set reference cell for field " << field.name()
                    << nl << "    Reference point " << refPointName
                    << " " << refPointi
                    << " found on " << sumHasRef << " domains (should be one)"
                    << nl << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unable to set reference cell for field " << field.name()
                << nl
                << "    Please supply either " << refCellName
                << " or " << refPointName << nl << exit(FatalIOError);
        }

        refValue = dict.lookup<scalar>(refValueName);

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    return setRefCell(field, field, dict, refCelli, refValue, forceReference);
}


Foam::scalar Foam::getRefCellValue
(
    const volScalarField& field,
    const label refCelli
)
{
    scalar refCellValue = (refCelli >= 0 ? field[refCelli] : 0.0);
    return returnReduce(refCellValue, sumOp<scalar>());
}


// ************************************************************************* //
