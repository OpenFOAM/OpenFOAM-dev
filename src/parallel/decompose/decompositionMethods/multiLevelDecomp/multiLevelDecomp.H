/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::multiLevelDecomp

Description
    Decomposition given using consecutive application of decomposers.

SourceFiles
    multiLevelDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef multiLevelDecomp_H
#define multiLevelDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multiLevelDecomp Declaration
\*---------------------------------------------------------------------------*/

class multiLevelDecomp
:
    public decompositionMethod
{
    // Private data

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

        //- Disallow default bitwise copy construct and assignment
        void operator=(const multiLevelDecomp&);
        multiLevelDecomp(const multiLevelDecomp&);


public:

    //- Runtime type information
    TypeName("multiLevel");


    // Constructors

        //- Construct given the decomposition dictionary
        multiLevelDecomp(const dictionary& decompositionDict);


    //- Destructor
    virtual ~multiLevelDecomp()
    {}


    // Member Functions

        //- Is method parallel aware (i.e. does it synchronize domains across
        //  proc boundaries)
        virtual bool parallelAware() const;

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
            const pointField& cc,
            const scalarField& cWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
