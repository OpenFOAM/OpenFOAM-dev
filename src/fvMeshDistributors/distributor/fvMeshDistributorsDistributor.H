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
    Foam::fvMeshDistributors::distributor

Description
    Dynamic mesh redistribution using the distributor specified in
    decomposeParDict

Usage
    Example of single field based refinement in all cells:
    \verbatim
    distributor
    {
        type            distributor;

        libs            ("libfvMeshDistributors.so");

        // How often to redistribute
        redistributionInterval  10;

        // Maximum fractional cell distribution imbalance
        // before rebalancing
        maxImbalance    0.1;
    }
    \endverbatim

SourceFiles
    fvMeshDistributorsDistributor.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDistributorsDistributor_H
#define fvMeshDistributorsDistributor_H

#include "fvMeshDistributor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class decompositionMethod;

namespace fvMeshDistributors
{

/*---------------------------------------------------------------------------*\
                          Class distributor Declaration
\*---------------------------------------------------------------------------*/

class distributor
:
    public fvMeshDistributor
{
protected:

    // Protected Member Data

        //- Cache the decomposer/distributor
        autoPtr<decompositionMethod> distributor_;

        //- Time-step interval between redistribution calls
        label redistributionInterval_;

        //- Maximum imbalance between the ideal number of cells per processor
        //  and the maximum or minimum as a ratio mag(1 - nCells/idealNcells)
        scalar maxImbalance_;

        //- The time index used for updating
        label timeIndex_;


    // Protected Member Functions

        //- Read the projection parameters from dictionary
        void readDict();

        //- Distribute the mesh and mesh data
        void distribute(const labelList& distribution);


public:

    //- Runtime type information
    TypeName("distributor");


    // Constructors

        //- Construct from fvMesh
        explicit distributor(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        distributor(const distributor&) = delete;


    //- Destructor
    virtual ~distributor();


    // Member Functions

        //- Distribute the
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Writing

        //- Write using given format, version and compression
        virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const distributor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshDistributors
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
