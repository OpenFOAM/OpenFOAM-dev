/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::cellSets

Description

\*---------------------------------------------------------------------------*/

#ifndef cellSets_H
#define cellSets_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cellSets Declaration
\*---------------------------------------------------------------------------*/

class cellSets
{
public:

        label nTets;
        label nPyrs;
        label nPrisms;
        label nHexesWedges;
        label nPolys;

        labelList tets;
        labelList pyrs;
        labelList prisms;
        labelList wedges;
        labelList hexes;
        labelList polys;


    // Constructors

        //- Construct given the number of cells
        cellSets(const label nCells)
        :
            nTets(0),
            nPyrs(0),
            nPrisms(0),
            nHexesWedges(0),
            nPolys(0),

            tets(nCells),
            pyrs(nCells),
            prisms(nCells),
            wedges(nCells),
            hexes(nCells),
            polys(nCells)
        {}


    // Member Functions

        void setSize(const label nCells)
        {
            nTets = 0;
            nPyrs = 0;
            nPrisms = 0;
            nHexesWedges = 0;
            nPolys = 0;

            tets.setSize(nCells);
            pyrs.setSize(nCells);
            prisms.setSize(nCells);
            wedges.setSize(nCells);
            hexes.setSize(nCells);
            polys.setSize(nCells);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
