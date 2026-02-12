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
    Foam::functionObjects::rigidBodySectionalForceProbes

Description
    This function calculates the fluid and acceleration forces and moments at a
    number of section-planes through a given rigid body. It writes the result
    into a log file.

    The cut-planes are defined as normal to an axis in the local coordinate
    system of the rigid body, and through a number of coordinates along that
    axis.

    Note that this function does not take into account the effect of joints or
    restraints attached to the body, so if these are applied then the forces
    generated will be somewhat incomplete. Note that ramping, acceleration
    relaxation, and solution tolerance can also result in components of the
    forces and moments to be missing from the output.

    Example of function object specification
    \verbatim
    rigidBodySectionalForceProbes
    {
        type            rigidBodySectionalForceProbes;
        libs            ("librigidBodyForces.so");
        body            beam;
        patches         (beam);
        axis            x;
        coordinates     (4.2 5 5.8);
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required   | Default value
        type         | Type name: rigidBodySectionalForceProbes | yes |
        body         | The name of the rigid body | yes |
        patches      | Patches included in the forces calculation | yes |
        axis         | Local axis normal to the section planes | yes |
        coordinates  | Coordinates of the section planes | yes |
        p            | Pressure field name     | no         | p
        U            | Velocity field name     | no         | U
        rho          | Density field name      | no         | rho
        phase        | Phase name for phase-fraction | no   |
        rhoRef       | Reference density       | if incompressible |
        pRef         | Reference pressure      | no         | 0
    \endtable

SourceFiles
    rigidBodySectionalForceProbes.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodySectionalForceProbes_functionObject_H
#define rigidBodySectionalForceProbes_functionObject_H

#include "rigidBodySectionalForcesBase.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rigidBodySectionalForceProbes Declaration
\*---------------------------------------------------------------------------*/

class rigidBodySectionalForceProbes
:
    public rigidBodySectionalForcesBase,
    public logFiles
{
    // Private Data

        //- Coordinates of the section planes
        List<scalar> coordinates_;

        //- The sorted order of the coordinates
        labelList coordinateOrders_;


    // Private Member Functions

        //- Return the distances from the origin to the section planes
        virtual tmp<scalarField> distances() const;

        //- Output file header information
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("rigidBodySectionalForceProbes");


    // Constructors

        //- Construct from Time and dictionary
        rigidBodySectionalForceProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodySectionalForceProbes
        (
            const rigidBodySectionalForceProbes&
        ) = delete;


    //- Destructor
    virtual ~rigidBodySectionalForceProbes();


    // Member Functions

        //- Read the rigidBodySectionalForceProbes data
        virtual bool read(const dictionary&);

        //- Write the rigidBodySectionalForceProbes
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodySectionalForceProbes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
