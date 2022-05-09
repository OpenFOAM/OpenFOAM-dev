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
    Foam::processorCyclicFvPatch

Description
    Processor patch.

SourceFiles
    processorCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorCyclicFvPatch_H
#define processorCyclicFvPatch_H

#include "processorCyclicPolyPatch.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class processorCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class processorCyclicFvPatch
:
    public processorFvPatch
{
    // Private Data

        const processorCyclicPolyPatch& procPolyPatch_;

public:

    //- Runtime type information
    TypeName(processorCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorCyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            processorFvPatch(patch, bm),
            procPolyPatch_(refCast<const processorCyclicPolyPatch>(patch))
        {}


    // Member Functions

        //- Return message tag used for sending
        virtual int tag() const
        {
            // Allocate from Pstream?
            return procPolyPatch_.tag();
        }

        //- Return the referring patch ID
        label referPatchID() const
        {
            return procPolyPatch_.referPatchID();
        }

        //- Return a reference to the referring patch
        const cyclicFvPatch& referPatch() const
        {
            const fvPatch& fvp = this->boundaryMesh()[referPatchID()];
            return refCast<const cyclicFvPatch>(fvp);
        }

        //- Return transformation between the coupled patches
        virtual const transformer& transform() const
        {
            return procPolyPatch_.transform();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
