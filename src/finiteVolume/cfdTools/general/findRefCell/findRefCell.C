/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::setRefCell
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
                refCelli = readLabel(dict.lookup(refCellName));

                if (refCelli < 0 || refCelli >= field.mesh().nCells())
                {
                    FatalIOErrorIn
                    (
                        "void Foam::setRefCell\n"
                         "(\n"
                         "    const volScalarField&,\n"
                         "    const volScalarField&,\n"
                         "    const dictionary&,\n"
                         "    label& scalar&,\n"
                         "    bool\n"
                         ")",
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

            // Note: find reference cell using facePlanes to avoid constructing
            //       face decomposition structure. Most likely the reference
            //       cell is an undistorted one so this should not be a
            //       problem.

            refCelli = field.mesh().findCell
            (
                refPointi,
                polyMesh::FACEPLANES
            );
            label hasRef = (refCelli >= 0 ? 1 : 0);
            label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
            if (sumHasRef != 1)
            {
                FatalIOErrorIn
                (
                    "void Foam::setRefCell\n"
                     "(\n"
                     "    const volScalarField&,\n"
                     "    const volScalarField&,\n"
                     "    const dictionary&,\n"
                     "    label& scalar&,\n"
                     "    bool\n"
                     ")",
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
            FatalIOErrorIn
            (
                "void Foam::setRefCell\n"
                 "(\n"
                 "    const volScalarField&,\n"
                 "    const volScalarField&,\n"
                 "    const dictionary&,\n"
                 "    label& scalar&,\n"
                 "    bool\n"
                 ")",
                dict
            )   << "Unable to set reference cell for field " << field.name()
                << nl
                << "    Please supply either " << refCellName
                << " or " << refPointName << nl << exit(FatalIOError);
        }

        refValue = readScalar(dict.lookup(refValueName));
    }
}


void Foam::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    setRefCell(field, field, dict, refCelli, refValue, forceReference);
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
