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
    Foam::randomDecomp

Description
    Random decomposition. Good for testing. Very bad for anything else.

SourceFiles
    randomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef randomDecomp_H
#define randomDecomp_H

#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class randomDecomp Declaration
\*---------------------------------------------------------------------------*/

class randomDecomp
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
        randomDecomp(const dictionary& decompositionDict);

        //- Disallow default bitwise copy construction
        randomDecomp(const randomDecomp&) = delete;


    //- Destructor
    virtual ~randomDecomp()
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
        void operator=(const randomDecomp&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
