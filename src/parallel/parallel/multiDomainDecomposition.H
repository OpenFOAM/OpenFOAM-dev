/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::multiDomainDecomposition

Description
    ...

SourceFiles
    multiDomainDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef multiDomainDecomposition_H
#define multiDomainDecomposition_H

#include "domainDecomposition.H"
#include "MultiRegionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class multiDomainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class multiDomainDecomposition
:
    public MultiRegionList<domainDecomposition>
{
    // Private Data

        //- Run times
        const processorRunTimes& runTimes_;


    // Private Member Functions

        //- Initialise and return the list of regions
        Foam::PtrList<domainDecomposition> init
        (
            const processorRunTimes& runTimes,
            const fileName& meshPath,
            const wordList& regionNames
        );


public:

    //- Runtime type information
    TypeName("multiDomainDecomposition");


    // Constructors

        //- Construct from processor run times and region names
        multiDomainDecomposition
        (
            const processorRunTimes& procRunTimes,
            const fileName& meshPath,
            const wordList& regionNames
        );


    //- Destructor
    virtual ~multiDomainDecomposition();


    // Member Functions

        //- Return the number of processors in the decomposition
        inline label nProcs() const
        {
            return runTimes_.nProcs();
        }

        //- Read in the complete mesh. Read the processor meshes if they are
        //  available and up to date relative to the complete mesh, or
        //  decompose if not. Return whether or not decomposition happened.
        bool readDecompose(const bool doSets);

        //- Read in the processor meshes. Read the complete mesh if it is
        //  available and up to date relative to the processor meshes, or
        //  reconstruct if not. Return whether or not reconstruction happened.
        bool readReconstruct(const bool doSets);

        //- Read-update for decomposition
        fvMesh::readUpdateState readUpdateDecompose();

        //- Read-update for reconstruction
        fvMesh::readUpdateState readUpdateReconstruct();

        //- Write the decomposed meshes and associated data
        void writeComplete(const bool doSets) const;

        //- Write the decomposed meshes and associated data
        void writeProcs(const bool doSets) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
