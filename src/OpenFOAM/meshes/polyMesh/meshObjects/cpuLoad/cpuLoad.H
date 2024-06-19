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
    Foam::cpuLoad

Description
    Class to maintain a field of the CPU load per cell.

    The optionalCpuLoad is a base-class for cpuLoad to simplify the
    implementation of optional CPU time caching.  This is achieved via the
    optionalCpuLoad::New function which returns a optionalCpuLoad with dummy
    functions if loadBalancing is false otherwise it creates or looks-up and
    returns a cpuLoad with the given name.

    Used for loadBalancing.

SourceFiles
    cpuLoad.C

\*---------------------------------------------------------------------------*/

#ifndef cpuLoad_H
#define cpuLoad_H

#include "cpuTime.H"
#include "scalarField.H"
#include "polyMesh.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class optionalCpuLoad Declaration
\*---------------------------------------------------------------------------*/

class optionalCpuLoad
{
    // Private Member Data

        //- Dummy optionalCpuLoad
        //  returned by optionalCpuLoad::New if loadBalancing = false
        static optionalCpuLoad optionalCpuLoad_;


public:

    // Constructors

        //- Construct from mesh, name and switch
        optionalCpuLoad()
        {}

        //- Disallow default bitwise copy construction
        optionalCpuLoad(const optionalCpuLoad&) = delete;


    // Selectors


        //- Construct from polyMesh if loadBalancing is true
        //  otherwise return the dummy optionalCpuLoad
        static optionalCpuLoad& New
        (
            const word& name,
            const polyMesh& mesh,
            const bool loadBalancing
        );


    //- Destructor
    virtual ~optionalCpuLoad()
    {}


    // Member Functions

        //- Reset the CPU time (dummy)
        virtual void resetCpuTime()
        {}

        //- Cache the CPU time increment for celli (dummy)
        virtual void cpuTimeIncrement(const label celli)
        {}

        //- Reset the CPU load field (dummy)
        virtual void reset()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const optionalCpuLoad&) = delete;
};


/*---------------------------------------------------------------------------*\
                          Class cpuLoad Declaration
\*---------------------------------------------------------------------------*/

class cpuLoad
:
    public optionalCpuLoad,
    public DemandDrivenMeshObject<polyMesh, TopoChangeableMeshObject, cpuLoad>,
    public scalarField
{
    // Private Data

        cpuTime cpuTime_;


public:

    // Declare name of the class and its debug switch
    TypeName("cpuLoad");


    // Constructors

        //- Construct from mesh, name and switch
        cpuLoad
        (
            const word& name,
            const polyMesh& mesh
        );

        //- Disallow default bitwise copy construction
        cpuLoad(const cpuLoad&) = delete;


    //- Destructor
    virtual ~cpuLoad();


    // Member Functions

        //- Reset the CPU time
        virtual void resetCpuTime();

        //- Cache the CPU time increment for celli
        virtual void cpuTimeIncrement(const label celli);

        //- Reset the CPU load field
        virtual void reset();


        //- Update for mesh motion
        virtual bool movePoints();

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap& map);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap& map);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cpuLoad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
