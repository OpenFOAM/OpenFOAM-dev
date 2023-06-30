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
    Foam::processorCyclicPointPatch

Description
    Processor patch boundary needs to be such that the ordering of
    points in the patch is the same on both sides.

    Looking at the creation of the faces on both sides of the processor
    patch they need to be identical on both sides with the normals pointing
    in opposite directions.  This is achieved by calling the reverseFace
    function in the decomposition.  It is therefore possible to re-create
    the ordering of patch points on the slave side by reversing all the
    patch faces of the owner.

SourceFiles
    processorCyclicPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicPointPatch_H
#define processorCyclicPointPatch_H

#include "processorPointPatch.H"
#include "processorCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class processorCyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicPointPatch
:
    public processorPointPatch
{
    // Private Data

        const processorCyclicPolyPatch& procCycPolyPatch_;


public:

    //- Runtime type information
    TypeName(processorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );

        //- Disallow default bitwise copy construction
        processorCyclicPointPatch(const processorCyclicPointPatch&) = delete;


    // Destructor

        virtual ~processorCyclicPointPatch();


    // Member Functions


        //- Return message tag to use for communication
        virtual int tag() const
        {
            return procCycPolyPatch_.tag();
        }

        //- Return the underlying processorCyclicPolyPatch
        const processorCyclicPolyPatch& procCyclicPolyPatch() const
        {
            return procCycPolyPatch_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const processorCyclicPointPatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
