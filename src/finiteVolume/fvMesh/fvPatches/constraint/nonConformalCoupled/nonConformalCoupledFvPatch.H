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
    Foam::nonConformalCoupledFvPatch

Description
    Non-conformal coupled FV patch. As nonConformalFvPatch, but is also coupled
    to another non-conformal patch.

See also
    Foam::nonConformalCoupledPolyPatch
    Foam::nonConformalFvPatch

SourceFiles
    nonConformalCoupledFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCoupledFvPatch_H
#define nonConformalCoupledFvPatch_H

#include "cyclicFvPatch.H"
#include "nonConformalCoupledPolyPatch.H"
#include "nonConformalFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonConformalErrorFvPatch;

/*---------------------------------------------------------------------------*\
                  Class nonConformalCoupledFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCoupledFvPatch
:
    public nonConformalFvPatch
{
    // Private Data

        //- Reference to the fvPatch
        const coupledFvPatch& patch_;

        //- Reference to the polyPatch
        const nonConformalCoupledPolyPatch& nonConformalCoupledPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights
        (
            scalarField& w,
            const vectorField& nbrSf,
            const vectorField& nbrDelta
        ) const;


public:

    //- Runtime type information
    TypeName(nonConformalCoupledPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch
        nonConformalCoupledFvPatch(const fvPatch& patch);


    //- Destructor
    virtual ~nonConformalCoupledFvPatch();


    // Member Functions

        // Access

            //- Return the reference to the polyPatch
            const nonConformalCoupledPolyPatch&
                nonConformalCoupledPatch() const;

            //- Does this side own the patch?
            bool owner() const;

            //- Does the coupled side own the patch?
            bool neighbour() const;

            //- Return transformation between the coupled patches
            const transformer& transform() const;

            //- Error patch name
            const word& errorPatchName() const;

            //- Error patch ID
            label errorPatchID() const;

            //- Error patch
            const nonConformalErrorFvPatch& errorPatch() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
