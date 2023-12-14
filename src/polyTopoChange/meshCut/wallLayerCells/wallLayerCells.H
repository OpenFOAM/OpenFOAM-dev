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

Class
    Foam::wallLayerCells

Description
    After construction is a list of \<refineCell\> of cells within a certain
    distance to the wall and with the wall normal direction.

SourceFiles
    wallLayerCells.C

\*---------------------------------------------------------------------------*/

#ifndef wallLayerCells_H
#define wallLayerCells_H

#include "edgeVertex.H"
#include "refineCell.H"
#include "typeInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class wallLayerCells Declaration
\*---------------------------------------------------------------------------*/

class wallLayerCells
:
    public edgeVertex,
    public List<refineCell>
{
    // Private Member Functions

        //- Check if any of the faces of celli is on processorPatch
        bool usesCoupledPatch(const label celli) const;

public:

    //- Runtime type information
    ClassName("wallLayerCells");


    // Constructors

        //- Construct from components
        wallLayerCells
        (
            const polyMesh& mesh,
            const List<word>& patchNames,
            const label nLayers
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
