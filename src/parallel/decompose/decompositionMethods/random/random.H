/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::decompositionMethods::random

Description
    Random decomposition. Good for testing. Very bad for anything else.

SourceFiles
    random.C

\*---------------------------------------------------------------------------*/

#ifndef random_H
#define random_H

#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                        Class random Declaration
\*---------------------------------------------------------------------------*/

class random
:
    public decompositionMethod
{
    // Private Data

        //- Generator seed
        label seed_;


public:

    //- Runtime type information
    TypeName("random");


    // Constructors

        //- Construct given the decomposition dictionary
        random(const dictionary& decompositionDict);

        //- Disallow default bitwise copy construction
        random(const random&) = delete;


    //- Destructor
    virtual ~random()
    {}


    // Member Functions

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity. Does not use the mesh.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        )
        {
            NotImplemented;
            return labelList(0);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const random&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
