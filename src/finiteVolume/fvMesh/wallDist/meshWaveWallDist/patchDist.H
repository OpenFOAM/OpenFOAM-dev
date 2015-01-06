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
    Foam::patchDist

Description
    Calculation of distance to nearest patch for all cells and boundary
    using meshWave.

    Distance correction (correctWalls = true):
        For each cell with face on wall calculate the true nearest point (by
        triangle decomposition) on that face and do the same for that face's
        pointNeighbours. This will find the true nearest distance in almost all
        cases. Only very skewed cells or cells close to another wall might be
        missed.

        For each cell with only one point on wall the same is done except now it
        takes the pointFaces() of the wall point to look for the nearest point.

SourceFiles
    patchDist.C

\*---------------------------------------------------------------------------*/

#ifndef patchDist_H
#define patchDist_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class patchDist Declaration
\*---------------------------------------------------------------------------*/

class patchDist
{
    // Private Member Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Set of patch IDs
        const labelHashSet patchIDs_;

        //- Do accurate distance calculation for near-wall cells.
        const bool correctWalls_;

        //- Number of unset cells and faces.
        mutable label nUnset_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchDist(const patchDist&);

        //- Disallow default bitwise assignment
        void operator=(const patchDist&);


public:

    // Constructors

        //- Construct from mesh and flag whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        patchDist
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs,
            const bool correctWalls = true
        );


    // Member Functions

        //- Calculate and return the distance to nearest patch
        //  for all cells and boundary
        tmp<volScalarField> y() const;

        label nUnset() const
        {
            return nUnset_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
