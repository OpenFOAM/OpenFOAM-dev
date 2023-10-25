/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "primitiveMeshCheck.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::meshCheck::checkPointNearness
(
    const primitiveMesh& mesh,
    const bool report,
    const scalar reportDistSqr,
    labelHashSet* setPtr
)
{
    const pointField& points = mesh.points();

    // Sort points
    SortableList<scalar> sortedMag(magSqr(points));

    label nClose = 0;

    for (label i = 1; i < sortedMag.size(); i++)
    {
        label pti = sortedMag.indices()[i];

        // Compare pti to any previous points with similar sortedMag
        for
        (
            label j = i-1;
            j >= 0 && (sortedMag[j] > sortedMag[i]-reportDistSqr);
            --j
        )
        {
            label prevPtI = sortedMag.indices()[j];

            if (magSqr(points[pti] - points[prevPtI]) < reportDistSqr)
            {
                nClose++;

                if (setPtr)
                {
                    setPtr->insert(pti);
                    setPtr->insert(prevPtI);
                }
            }
        }
    }

    reduce(nClose, sumOp<label>());

    if (nClose > 0)
    {
        if (report)
        {
            Info<< "  <<Points closer than " << Foam::sqrt(reportDistSqr)
                << " together found, number: " << nClose
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
