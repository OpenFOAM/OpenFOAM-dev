/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::nonConformalErrorPolyPatch

Description
    Non-conformal error poly patch. As nonConformalPolyPatch. This patch is a
    non-coupled non-conformal patch which is used to manage the errors created
    during the construction of a non-conformal coupled interface. Every patch
    used as the original patch of the owner side of a non-conformal coupled
    interface must also have an associated error patch.

See also
    Foam::nonConformalPolyPatch

SourceFiles
    nonConformalErrorPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorPolyPatch_H
#define nonConformalErrorPolyPatch_H

#include "polyPatch.H"
#include "nonConformalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalCyclicPolyPatch;

/*---------------------------------------------------------------------------*\
                 Class nonConformalErrorPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalErrorPolyPatch
:
    public polyPatch,
    public nonConformalPolyPatch
{
protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex);


public:

    //- Runtime type information
    TypeName("nonConformalError");


    // Constructors

        //- Construct from components
        nonConformalErrorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from components
        nonConformalErrorPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& origPatchName
        );

        //- Construct from dictionary
        nonConformalErrorPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        nonConformalErrorPolyPatch
        (
            const nonConformalErrorPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        nonConformalErrorPolyPatch
        (
            const nonConformalErrorPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& origPatchName
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new nonConformalErrorPolyPatch(*this, bm)
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
                new nonConformalErrorPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    origPatchName()
                )
            );
        }


    //- Destructor
    virtual ~nonConformalErrorPolyPatch();


    // Member Functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
