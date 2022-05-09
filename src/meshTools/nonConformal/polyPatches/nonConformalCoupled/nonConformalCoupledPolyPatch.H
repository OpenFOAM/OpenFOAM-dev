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
    Foam::nonConformalCoupledPolyPatch

Description
    Non-conformal coupled poly patch. As nonConformalPolyPatch, but this patch
    is coupled to another non-conformal patch. Information about the
    owner/neighbour relationship and transformation are made available by this
    class. Also, this patch allows access to the error patch.

See also
    Foam::nonConformalPolyPatch
    Foam::nonConformalErrorPolyPatch

SourceFiles
    nonConformalCoupledPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCoupledPolyPatch_H
#define nonConformalCoupledPolyPatch_H

#include "cyclicPolyPatch.H"
#include "nonConformalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalErrorPolyPatch;

/*---------------------------------------------------------------------------*\
                 Class nonConformalCoupledPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCoupledPolyPatch
:
    public nonConformalPolyPatch
{
private:

    // Private

        //- Reference to the coupled polyPatch
        const coupledPolyPatch& patch_;

        //- Name of the error patch associated with this cyclic
        mutable word errorPatchName_;

        //- Index of the error patch associated with this cyclic
        mutable label errorPatchID_;


protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex);


public:

    //- Runtime type information
    TypeName("nonConformalCoupled");


    // Constructors

        //- Construct from a patch
        nonConformalCoupledPolyPatch(const polyPatch& patch);

        //- Construct from a patch and an original patch name
        nonConformalCoupledPolyPatch
        (
            const polyPatch& patch,
            const word& origPatchName
        );

        //- Construct from a patch and a dictionary
        nonConformalCoupledPolyPatch
        (
            const polyPatch& patch,
            const dictionary& dict
        );

        //- Construct from a non-conformal coupled patch
        nonConformalCoupledPolyPatch(const nonConformalCoupledPolyPatch& patch);


    //- Destructor
    virtual ~nonConformalCoupledPolyPatch();


    // Member Functions

        //- Does this side own the patch?
        bool owner() const;

        //- Does the coupled side own the patch?
        bool neighbour() const;

        //- Return transformation between the coupled patches
        const transformer& transform() const;

        //- Error patch name
        const word& errorPatchName() const;

        //- Error patchID
        label errorPatchID() const;

        //- Error patch
        const nonConformalErrorPolyPatch& errorPatch() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
