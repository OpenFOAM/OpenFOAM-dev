/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::rigidBodyForces

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches of a moving rigid body.

    The centre of rotation (CofR) of the moving rigid object is obtained
    directly from the corresponding Foam::RBD::rigidBodyMotion of the
    specified body.

    Member function rigidBodyForces::write() calculates the forces/moments and
    writes the forces/moments into the file \<timeDir\>/rigidBodyForces.dat
    and bin data (if selected) to the file \<timeDir\>/rigidBodyForces_bin.dat

    Example of function object specification:
    \verbatim
    rigidBodyForces1
    {
        type        rigidBodyForces;
        libs        ("librigidBodyForces.so");

        body        (hull);
        patches     (walls);

        log         yes;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: rigidBodyForces       | yes         |
        log          | Write force data to standard output | no | no
        body         | Name of the rigid body  | yes |
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        phase        | Phase name for phase-fraction  | no   |
        directForceDensity | Force density supplied directly (see below)|no|no
        fD           | Name of force density field (see below) | no | fD
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according o
    \table
        nBin         | number of data bins     | yes         |
        direction    | direction along which bins are defined | yes |
        cumulative   | bin data accumulated with increasing distance | yes |
    \endtable

  Note
    - For incompressible cases, set \c rho to \c rhoInf and provide
      a \c rhoInf value corresponding to the free-stream constant density.
    - If the \c phase name is specified the corresponding phase-fraction field
      \c alpha.<phase> is used to filter the surface force field
      before integration.
    - If the force density is supplied directly, set the \c directForceDensity
      flag to 'yes', and supply the force density field using the \c
      fDName entry

See also
    Foam::functionObject
    Foam::functionObjects::forcesBase
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::timeControl
    Foam::RBD::rigidBodyMotion

SourceFiles
    rigidBodyForces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_rigidBodyForces_H
#define functionObjects_rigidBodyForces_H

#include "forcesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class rigidBodyForces Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyForces
:
    public forcesBase
{
    // Private Data

        //- Name of the rigid body
        word body_;


protected:

    // Protected Member Functions

        //- Return the current centre of the rigid body
        virtual vector CofR() const;


public:

    //- Runtime type information
    TypeName("rigidBodyForces");


    // Constructors

        //- Construct from Time and dictionary
        rigidBodyForces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        rigidBodyForces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        rigidBodyForces(const rigidBodyForces&) = delete;


    //- Destructor
    virtual ~rigidBodyForces();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyForces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
