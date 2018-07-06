/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::patchDistMethods::meshWave

Description
    Fast topological mesh-wave method for calculating the distance to nearest
    patch for all cells and boundary.

    For regular/un-distorted meshes this method is accurate but for skewed,
    non-orthogonal meshes it is approximate with the error increasing with the
    degree of mesh distortion.  The distance from the near-wall cells to the
    boundary may optionally be corrected for mesh distortion by setting
    correctWalls = true.

    Example of the wallDist specification in fvSchemes:
    \verbatim
        wallDist
        {
            method meshWave;

            // Optional entry enabling the calculation
            // of the normal-to-wall field
            nRequired false;
        }
    \endverbatim

See also
    Foam::patchDistMethod::Poisson
    Foam::wallDist

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

        //- Construct from coefficients dictionary, mesh
        //  and fixed-value patch set
        meshWave
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelHashSet& patchIDs
        );

        //- Construct from mesh, fixed-value patch set and flag specifying
        //  whether or not to correct wall.
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

        //- Correct the given distance-to-patch and normal-to-patch fields
        virtual bool correct(volScalarField& y, volVectorField& n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
