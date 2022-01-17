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
    Foam::fvMeshDistributors::loadBalancer

Description
    Dynamic mesh redistribution using the distributor specified in
    decomposeParDict

Usage
    Example of single field based refinement in all cells:
    \verbatim
    distributor
    {
        type            loadBalancer;

        libs            ("libfvMeshDistributors.so");

        // How often to redistribute
        redistributionInterval  10;

        // Maximum fractional cell distribution imbalance
        // before rebalancing
        maxImbalance    0.1;
    }
    \endverbatim

SourceFiles
    fvMeshDistributorsloadBalancer.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDistributorsLoadBalancer_H
#define fvMeshDistributorsLoadBalancer_H

#include "fvMeshDistributorsDistributor.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{

/*---------------------------------------------------------------------------*\
                          Class loadBalancer Declaration
\*---------------------------------------------------------------------------*/

class loadBalancer
:
    public distributor
{
    // Private Member Data

        //- CPU time consumed during the time-step
        cpuTime cpuTime_;

        //- Enable multi-constraint load-balancing in which separate weights
        //  are provided to the distrubutor for each of the CPU loads.
        //  When disabled the CPU loads are summed and a single weight per cell
        //  is provided to the distrubutor.
        //  Defaults to true.
        Switch multiConstraint_;


    // Private Member Functions

        //- Read the projection parameters from dictionary
        void readDict();


public:

    //- Runtime type information
    TypeName("loadBalancer");


    // Constructors

        //- Construct from fvMesh
        explicit loadBalancer(fvMesh& mesh);


    //- Destructor
    virtual ~loadBalancer();


    // Member Functions

        //- Distribute the
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshDistributors
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
