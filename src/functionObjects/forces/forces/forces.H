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
    Foam::functionObjects::forces

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches.

    Member function forces::write() calculates the forces/moments and
    writes the forces/moments into the file \<timeDir\>/forces.dat and bin
    data (if selected) to the file \<timeDir\>/forces_bin.dat

    Example of function object specification:
    \verbatim
    forces1
    {
        type        forces;
        libs        ("libforces.so");
        ...
        log         yes;
        patches     (walls);

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: forces       | yes         |
        log          | Write force data to standard output | no | no
        patches      | Patches included in the forces calculation | yes |
        p            | Pressure field name     | no          | p
        U            | Velocity field name     | no          | U
        rho          | Density field name (see below) | no   | rho
        phase        | Phase name for phase-fraction  | no   |
        CofR         | Centre of rotation      | yes  |
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
    Foam::forceCoeffs

SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forces_H
#define functionObjects_forces_H

#include "forcesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class forces
:
    public forcesBase
{

protected:

    // Protected data

            //- Centre of rotation
            vector CofR_;


    // Protected Member Functions

        //- Return the current centre of the rigid body
        virtual vector CofR() const;

        //- Write the constant centre of rotation value in the header
        virtual void writeCoRValueHeader(Ostream& file);

        //- Do not write the constant centre of rotation
        virtual void writeCoRHeader(Ostream& file);

        //- Do not write the constant centre of rotation
        virtual void writeCofR(Ostream& file);


public:

    //- Runtime type information
    TypeName("forces");


    // Constructors

        //- Construct from Time and dictionary
        forces
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        forces
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        forces(const forces&) = delete;


    //- Destructor
    virtual ~forces();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const forces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
