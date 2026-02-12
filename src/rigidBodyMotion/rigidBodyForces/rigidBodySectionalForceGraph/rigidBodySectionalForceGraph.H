/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::functionObjects::rigidBodySectionalForceGraph

 Description
    This function calculates the fluid and acceleration forces and moments at a
    series of section-planes through a given rigid body. It writes the result
    into a graph which can be used to construct (potentially just a part of) a
    shear force and bending moment diagram. Points in the graph are graded to
    approximately match the resolution of the mesh.

    The cut-planes are defined as normal to an axis in the local coordinate
    system of the rigid body, and the line of the graph extends along the axis
    from an origin defined in the same local coordinate system. Forces and
    moments are calculated on the side of the plane on which coordinates have a
    positive axis-component.

    Note that this function does not take into account the effect of joints or
    restraints attached to the body, so if these are applied then the shear
    force and bending moment diagrams generated will be somewhat incomplete.
    Note that ramping, acceleration relaxation, and solution tolerance can also
    result in components of the forces and moments to be missing from the
    output.

    Example of function object specification:
    \verbatim
    rigidBodySectionalForceGraph
    {
        type            rigidBodySectionalForceGraph;
        libs            ("librigidBodyForces.so");
        body            beam;
        patches         (beam);
        axis            x;
        nPoints         101;
        setFormat       raw;
        writeControl    writeTime;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required   | Default value
        type         | Type name: rigidBodySectionalForceGraph | yes |
        body         | The name of the rigid body | yes |
        patches      | Patches included in the forces calculation | yes |
        axis         | Local axis normal to the section planes | yes |
        origin       | Local origin of the plot | no        | body centroid
        nPoints      | Number of points in the plot | yes   |
        setFormat    | Graph output format     | yes        |
        p            | Pressure field name     | no         | p
        U            | Velocity field name     | no         | U
        rho          | Density field name      | no         | rho
        phase        | Phase name for phase-fraction | no   |
        rhoRef       | Reference density       | if incompressible |
        pRef         | Reference pressure      | no         | 0
    \endtable

SourceFiles
    rigidBodySectionalForceGraph.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodySectionalForceGraph_functionObject_H
#define rigidBodySectionalForceGraph_functionObject_H

#include "rigidBodySectionalForcesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rigidBodySectionalForceGraph Declaration
\*---------------------------------------------------------------------------*/

class rigidBodySectionalForceGraph
:
    public rigidBodySectionalForcesBase
{
    // Private Data

        //- The distances from the origin to the section planes
        mutable autoPtr<scalarField> distancesPtr_;

        //- The number of points in the plot
        label nPoints_;

        //- File writer
        autoPtr<setWriter> formatter_;

        //- Number of optimisation iterations. Default is 2.
        label nOptimiseIter_;


    // Private Member Functions

        //- Clear the cached weights
        virtual void clear();

        //- Return the distances from the origin to the section planes
        virtual tmp<scalarField> distances() const;


public:

    //- Runtime type information
    TypeName("rigidBodySectionalForceGraph");


    // Constructors

        //- Construct from Time and dictionary
        rigidBodySectionalForceGraph
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodySectionalForceGraph
        (
            const rigidBodySectionalForceGraph&
        ) = delete;


    //- Destructor
    virtual ~rigidBodySectionalForceGraph();


    // Member Functions

        //- Read the rigidBodySectionalForceGraph data
        virtual bool read(const dictionary&);

        //- Write the rigidBodySectionalForceGraph
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodySectionalForceGraph&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
