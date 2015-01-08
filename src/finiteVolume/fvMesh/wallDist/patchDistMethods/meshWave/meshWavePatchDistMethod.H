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
    meshWavePatchDistMethod.C

\*---------------------------------------------------------------------------*/

#ifndef meshWavePatchDistMethod_H
#define meshWavePatchDistMethod_H

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace patchDistMethods
{

/*---------------------------------------------------------------------------*\
                          Class meshWave Declaration
\*---------------------------------------------------------------------------*/

class meshWave
:
    public patchDistMethod
{
    // Private Member Data

        //- Do accurate distance calculation for near-wall cells.
        const bool correctWalls_;

        //- Number of unset cells and faces.
        mutable label nUnset_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshWave(const meshWave&);

        //- Disallow default bitwise assignment
        void operator=(const meshWave&);


public:

    //- Runtime type information
    TypeName("meshWave");


    // Constructors

        //- Construct from mesh and flag whether or not to correct wall.
        //  Calculate for all cells.
        meshWave
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Construct from mesh and flag whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        meshWave
        (
            const fvMesh& mesh,
            const labelHashSet& patchIDs,
            const bool correctWalls = true
        );


    // Member Functions

        label nUnset() const
        {
            return nUnset_;
        }

        //- Correct the given distance-to-patch field
        virtual bool correct(volScalarField& y);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
