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
    Foam::nonConformalProcessorCyclicPolyPatch

Description
    Non-conformal processor cyclic poly patch. As nonConformalCyclicPolyPatch,
    but the neighbouring patch is on a different processor.

See also
    Foam::nonConformalCyclicPatch

SourceFiles
    nonConformalProcessorCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalProcessorCyclicPolyPatch_H
#define nonConformalProcessorCyclicPolyPatch_H

#include "processorCyclicPolyPatch.H"
#include "nonConformalCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class nonConformalProcessorCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalProcessorCyclicPolyPatch
:
    public processorCyclicPolyPatch,
    public nonConformalCoupledPolyPatch
{
protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex);


public:

    //- Runtime type information
    TypeName("nonConformalProcessorCyclic");


    // Constructors

        //- Construct from components
        nonConformalProcessorCyclicPolyPatch
        (
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const int myProcNo,
            const int neighbProcNo,
            const word& referPatchName,
            const word& origPatchName,
            const word& patchType = typeName
        );

        //- Construct from dictionary
        nonConformalProcessorCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        nonConformalProcessorCyclicPolyPatch
        (
            const nonConformalProcessorCyclicPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        nonConformalProcessorCyclicPolyPatch
        (
            const nonConformalProcessorCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new nonConformalProcessorCyclicPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new nonConformalProcessorCyclicPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor
    virtual ~nonConformalProcessorCyclicPolyPatch();


    // Member Functions

        //- Get the referred patch
        const nonConformalCyclicPolyPatch& referPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[referPatchID()];
            return refCast<const nonConformalCyclicPolyPatch>(pp);
        }

        //- Inherit the processor cyclic owner method
        using processorCyclicPolyPatch::owner;

        //- Inherit the processor cyclic neighbour method
        using processorCyclicPolyPatch::neighbour;

        //- Inherit the processor cyclic transform method
        using processorCyclicPolyPatch::transform;

        //- Is this patch coupled? Returns false. For NCC patches the poly
        //  mesh is considered non-coupled whilst the finite volume mesh is
        //  considered coupled.
        virtual bool coupled() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
