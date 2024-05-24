/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::decompositionMethods::multiLevel

Description
    Decomposition given using consecutive application of decomposers.

SourceFiles
    multiLevel.C

\*---------------------------------------------------------------------------*/

#ifndef multiLevel_H
#define multiLevel_H

#include "decompositionMethod.H"

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                      Class multiLevel Declaration
\*---------------------------------------------------------------------------*/

class multiLevel
:
    public decompositionMethod
{
    // Private Data

        dictionary methodsDict_;

        PtrList<decompositionMethod> methods_;


    // Private Member Functions

        //- Given connectivity across processors work out connectivity
        //  for a (consistent) subset
        void subsetGlobalCellCells
        (
            const label nDomains,
            const label domainI,
            const labelList& dist,

            const labelListList& cellCells,
            const labelList& set,
            labelListList& subCellCells,
            labelList& cutConnections
        ) const;

        //- Decompose level methodI without addressing
        void decompose
        (
            const labelListList& pointPoints,
            const pointField& points,
            const scalarField& pointWeights,
            const labelList& pointMap,  // map back to original points
            const label levelI,

            labelField& finalDecomp
        );


public:

    //- Runtime type information
    TypeName("multiLevel");


    // Constructors

        //- Construct given the decomposition dictionary
        multiLevel(const dictionary& decompositionDict);

        //- Disallow default bitwise copy construction
        multiLevel(const multiLevel&) = delete;


    //- Destructor
    virtual ~multiLevel()
    {}


    // Member Functions

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multiLevel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
