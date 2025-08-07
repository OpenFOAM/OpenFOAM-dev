/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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
    Foam::pointDist

Description
    Calculates the distance to the specified sets of patch and pointZone points
    or for all points

SourceFiles
    pointDist.C

\*---------------------------------------------------------------------------*/

#ifndef pointDist_H
#define pointDist_H

#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class pointDist Declaration
\*---------------------------------------------------------------------------*/

class pointDist
:
    public pointScalarField
{
    // Private Member Data

        //- Reference to pointField
        const pointField& points_;

        //- Set of start patch IDs for which the distance is 0
        const labelHashSet startPatchIDs_;

        //- Set of start pointZone IDs for which the distance is 0
        const labelHashSet startZoneIDs_;

        //- Set of end patch IDs for which the distance is great
        const labelHashSet endPatchIDs_;

        //- Set of end pointZone IDs for which the distance is great
        const labelHashSet endZoneIDs_;

        //- Maximum distance at which the point-edge wave is stopped
        const scalar maxDist_;


public:

    // Constructors

        //- Construct from mesh and set of start patches
        //  for which the distance is 0 and points
        pointDist
        (
            const pointMesh& pMesh,
            const labelHashSet& startPatchIDs,
            const pointField& points,
            const scalar maxDist = rootVGreat
        );

        //- Construct from mesh, set of start patches and pointZones
        //  for which the distance is 0
        //  and points
        pointDist
        (
            const pointMesh& pMesh,
            const labelHashSet& startPatchIDs,
            const labelHashSet& startZoneIDs,
            const labelHashSet& endPatchIDs,
            const labelHashSet& endZoneIDs,
            const pointField& points,
            const scalar maxDist = rootVGreat
        );

        //- Disallow default bitwise copy construction
        pointDist(const pointDist&) = delete;


    //- Destructor
    virtual ~pointDist();


    // Member Functions

        const pointScalarField& y() const
        {
            return *this;
        }

        //- Correct for mesh geom/topo changes
        void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointDist&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
