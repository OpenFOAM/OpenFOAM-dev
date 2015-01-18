/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "indexedCellOps.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<typename CellType>
Foam::label CGAL::indexedCellOps::dualVertexMasterProc(const CellType& c)
{
    if (!c->parallelDualVertex())
    {
        return -1;
    }

    // The master processor is the lowest numbered of the four on this tet.

    int masterProc = Foam::Pstream::nProcs() + 1;

    for (Foam::label vI = 0; vI < 4; ++vI)
    {
        if (c->vertex(vI)->referred())
        {
            masterProc = min(masterProc, c->vertex(vI)->procIndex());
        }
        else
        {
            masterProc = min(masterProc, Foam::Pstream::myProcNo());
        }
    }

    return masterProc;
}


template<typename CellType>
Foam::FixedList<Foam::label, 4>
CGAL::indexedCellOps::processorsAttached(const CellType& c)
{
    Foam::FixedList<Foam::label, 4> procsAttached(Foam::Pstream::myProcNo());

    if (!c->parallelDualVertex())
    {
        return procsAttached;
    }

    for (Foam::label vI = 0; vI < 4; ++vI)
    {
        if (c->vertex(vI)->referred())
        {
            procsAttached[vI] = c->vertex(vI)->procIndex();
        }
    }

    return procsAttached;
}


// ************************************************************************* //
