/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::meshWaveWallDist

Description
    Specialisation of patchDist for wall distance calculation

SourceFiles
    meshWaveWallDist.C

\*---------------------------------------------------------------------------*/

#ifndef meshWaveWallDist_H
#define meshWaveWallDist_H

#include "patchDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class meshWaveWallDist Declaration
\*---------------------------------------------------------------------------*/

class meshWaveWallDist
:
    public patchDist
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshWaveWallDist(const meshWaveWallDist&);

        //- Disallow default bitwise assignment
        void operator=(const meshWaveWallDist&);


public:

    // Constructors

        //- Construct from mesh and flag whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        meshWaveWallDist
        (
            const fvMesh& mesh,
            const bool correctWalls = true
        );


    //- Destructor
    virtual ~meshWaveWallDist();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
