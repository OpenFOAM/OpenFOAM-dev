/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::decompositionMethods::none

Description
    Dummy decomposition method

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef none_H
#define none_H

#include "decompositionMethod.H"

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                          Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public decompositionMethod
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct given the decomposition dictionary
        none(const dictionary& decompositionDict);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    virtual ~none()
    {}


    // Member Functions

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cellCentres,
            const scalarField& cellWeights
        )
        {
            return labelList(cellCentres.size(), Pstream::myProcNo());
        }

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        )
        {
            return labelList(globalCellCells.size(), Pstream::myProcNo());
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
