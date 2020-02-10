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
    Foam::cyclicPointPatch

Description
    Cyclic patch for post-processing.

SourceFiles
    cyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPointPatch_H
#define cyclicPointPatch_H

#include "coupledFacePointPatch.H"
#include "cyclicPolyPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPointPatch
:
    public coupledFacePointPatch
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const cyclicPolyPatch& cyclicPolyPatch_;


    // Private Member Functions

        //- Disallow default bitwise copy construction
        cyclicPointPatch(const cyclicPointPatch&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicPointPatch&);


    // Demand driven private data

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName(cyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~cyclicPointPatch();


    // Member Functions

        // Access

            //- Return the constraint type this pointPatch implements.
            virtual const word& constraintType() const
            {
                return type();
            }

            //- Return the underlying cyclicPolyPatch
            const cyclicPolyPatch& cyclicPatch() const
            {
                return cyclicPolyPatch_;
            }

            //- Return neighbour point patch
            const cyclicPointPatch& nbrPatch() const
            {
                label patchi = cyclicPolyPatch_.nbrPatchID();
                const pointPatch& pp = this->boundaryMesh()[patchi];
                return refCast<const cyclicPointPatch>(pp);
            }

            //- Return transformation between the coupled patches
            virtual const transformer& transform() const
            {
                return cyclicPolyPatch_.transform();
            }


        // Access functions for demand driven data

            //- Return the set of pairs of points that require transformation
            //  and/or mapping. First index is on this patch, second on the
            //  neighbour patch.
            virtual const edgeList& transformPairs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
