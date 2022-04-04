/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::cyclicAMIPointPatch

Description
    Cyclic AMI point patch - place holder only

SourceFiles
    cyclicAMIPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMIPointPatch_H
#define cyclicAMIPointPatch_H

#include "coupledFacePointPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cyclicAMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMIPointPatch
:
    public coupledFacePointPatch
{
    // Private Data

        //- Local reference cast into the cyclic AMI patch
        const cyclicAMIPolyPatch& cyclicAMIPolyPatch_;


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initCalcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initTopoChange(PstreamBuffers&);

        //- Update of the patch topology
        virtual void topoChange(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName(cyclicAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );

        //- Disallow default bitwise copy construction
        cyclicAMIPointPatch(const cyclicAMIPointPatch&);


    //- Destructor
    virtual ~cyclicAMIPointPatch();


    // Member Functions

        //- Return true if this patch field is coupled
        virtual bool coupled() const;

        //- Return the underlying cyclicAMIPolyPatch
        const cyclicAMIPolyPatch& cyclicAMIPatch() const
        {
            return cyclicAMIPolyPatch_;
        }

        //- Return neighbour point patch
        const cyclicAMIPointPatch& nbrPatch() const
        {
            label patchi = cyclicAMIPolyPatch_.nbrPatchID();
            const pointPatch& pp = this->boundaryMesh()[patchi];
            return refCast<const cyclicAMIPointPatch>(pp);
        }

        //- Return transformation between the coupled patches
        virtual const transformer& transform() const
        {
            return cyclicAMIPolyPatch_.transform();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cyclicAMIPointPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
