/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "extendedCentredCellToCellStencil.H"
#include "distributionMap.H"
#include "cellToCellStencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCentredCellToCellStencil::extendedCentredCellToCellStencil
(
    const cellToCellStencil& stencil
)
:
    extendedCellToCellStencil(stencil.mesh()),
    stencil_(stencil)
{
    // Calculate distribute map (also renumbers elements in stencil)
    List<Map<label>> compactMap(Pstream::nProcs());
    mapPtr_.reset
    (
        new distributionMap
        (
            stencil.globalNumbering(),
            stencil_,
            compactMap
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::extendedCentredCellToCellStencil::compact()
{
    boolList isInStencil(map().constructSize(), false);

    forAll(stencil_, celli)
    {
        const labelList& stencilCells = stencil_[celli];

        forAll(stencilCells, i)
        {
            isInStencil[stencilCells[i]] = true;
        }
    }

    mapPtr_().compact(isInStencil, Pstream::msgType());
}


// ************************************************************************* //
