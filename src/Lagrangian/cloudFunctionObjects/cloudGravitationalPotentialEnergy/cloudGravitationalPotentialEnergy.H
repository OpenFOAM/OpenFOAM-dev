/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::cloudGravitationalPotentialEnergy

Description
    Function to compute and write the gravitational potential energy field for
    a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudGravitationalPotentialEnergy1
    {
        type            cloudGravitationalPotentialEnergy;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudGravitationalPotentialEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudGravitationalPotentialEnergy_H
#define functionObjects_cloudGravitationalPotentialEnergy_H

#include "cloudLagrangianMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class cloudGravitationalPotentialEnergy Declaration
\*---------------------------------------------------------------------------*/

class cloudGravitationalPotentialEnergy
:
    public cloudLagrangianMeshFunctionObject
{
private:

    // Private Member Data

        //- Gravitational acceleration
        const dimensionedVector& g;


public:

    //- Runtime type information
    TypeName("cloudGravitationalPotentialEnergy");


    // Public Static Data

        //- The name of the gravitational potential energy field
        static const word GPEName_;


    // Constructors

        //- Construct from Time and dictionary
        cloudGravitationalPotentialEnergy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudGravitationalPotentialEnergy
        (
            const cloudGravitationalPotentialEnergy&
        ) = delete;


    //- Destructor
    virtual ~cloudGravitationalPotentialEnergy();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Hook before solution steps
        virtual void preSolve();

        //- Calculate and store the gravitational potential energy
        virtual bool execute();

        //- Write the gravitational potential energy
        virtual bool write();

        //- Clear the gravitational potential energy
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudGravitationalPotentialEnergy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
