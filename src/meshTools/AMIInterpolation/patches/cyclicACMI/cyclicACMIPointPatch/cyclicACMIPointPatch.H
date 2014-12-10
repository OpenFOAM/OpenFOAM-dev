/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::cyclicACMIPointPatch

Description
    Cyclic AMI point patch - place holder only

SourceFiles
    cyclicACMIPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicACMIPointPatch_H
#define cyclicACMIPointPatch_H

#include "coupledFacePointPatch.H"
#include "cyclicACMIPolyPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class cyclicACMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicACMIPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic AMI patch
        const cyclicACMIPolyPatch& cyclicACMIPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        cyclicACMIPointPatch(const cyclicACMIPointPatch&);

        //- Disallow default assignment
        void operator=(const cyclicACMIPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

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
    TypeName(cyclicACMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicACMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~cyclicACMIPointPatch();


    // Member Functions

        //- Is patch 'coupled'. Note that on AMI the geometry is not
        //  coupled but the fields are!
        virtual bool coupled() const
        {
            return false;
        }

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Return the underlying cyclicAMIPolyPatch
        const cyclicACMIPolyPatch& cyclicACMIPatch() const
        {
            return cyclicACMIPolyPatch_;
        }

        //- Return neighbour point patch
        const cyclicACMIPointPatch& neighbPatch() const
        {
            label patchI = cyclicACMIPolyPatch_.neighbPatchID();
            const pointPatch& pp = this->boundaryMesh()[patchI];
            return refCast<const cyclicACMIPointPatch>(pp);
        }

        //- Are the cyclic planes parallel
        bool parallel() const
        {
            return cyclicACMIPolyPatch_.parallel();
        }

        //- Return face transformation tensor
        const tensorField& forwardT() const
        {
            return cyclicACMIPolyPatch_.forwardT();
        }

        //- Return neighbour-cell transformation tensor
        const tensorField& reverseT() const
        {
            return cyclicACMIPolyPatch_.reverseT();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
