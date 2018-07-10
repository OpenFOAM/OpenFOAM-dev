/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::springRenumber

Description
    Use spring analogy - attract neighbouring cells according to the distance
    of their cell indices.

    // Maximum jump of cell indices. Is fraction of number of cells
    maxCo 0.1;

    // Limit the amount of movement; the fraction maxCo gets decreased
    // with every iteration.
    freezeFraction 0.9;

    // Maximum number of iterations
    maxIter 1000;

SourceFiles
    springRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef springRenumber_H
#define springRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class springRenumber Declaration
\*---------------------------------------------------------------------------*/

class springRenumber
:
    public renumberMethod
{
    // Private data

        const dictionary& dict_;

        const scalar maxCo_;

        const label maxIter_;

        const scalar freezeFraction_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const springRenumber&);
        springRenumber(const springRenumber&);


public:

    //- Runtime type information
    TypeName("spring");


    // Constructors

        //- Construct given the renumber dictionary
        springRenumber(const dictionary& renumberDict);


    //- Destructor
    virtual ~springRenumber()
    {}


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
            return labelList(0);
        }

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
