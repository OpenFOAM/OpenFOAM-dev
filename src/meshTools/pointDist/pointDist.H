/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

        //- Set of patch IDs
        const labelHashSet patchIDs_;

        //- Set of pointZone IDs for which the distance is 0
        const labelHashSet zoneIDs_;

        //- Number of unset points
        label nUnset_;


public:

    // Constructors

        //- Construct from mesh and set of patches and points
        pointDist
        (
            const pointMesh& pMesh,
            const labelHashSet& patchIDs,
            const pointField& points
        );

        //- Construct from mesh, set of patches,
        //  set of pointZones for which the distance is 0
        //  and points
        pointDist
        (
            const pointMesh& pMesh,
            const labelHashSet& patchIDs,
            const labelHashSet& zoneIDs,
            const pointField& points
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

        label nUnset() const
        {
            return nUnset_;
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
