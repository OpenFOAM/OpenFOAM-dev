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
    Foam::functionObjects::cloudBoundaryCollisionNumberFlux

Description
    Function to write the boundary collision number flux field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudBoundaryCollisionNumberFlux1
    {
        type            cloudBoundaryCollisionNumberFlux;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudBoundaryCollisionNumberFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudBoundaryCollisionNumberFlux_H
#define functionObjects_cloudBoundaryCollisionNumberFlux_H

#include "cloudBoundaryCollisionFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
               Class cloudBoundaryCollisionNumberFlux Declaration
\*---------------------------------------------------------------------------*/

class cloudBoundaryCollisionNumberFlux
:
    public cloudBoundaryCollisionFlux
{
private:

    // Private Member Functions

        //- Return the property for which to compute the flux
        virtual tmp<LagrangianSubScalarSubField> q
        (
            const LagrangianSubScalarSubField& fraction,
            const label sign
        ) const;


public:

    //- Runtime type information
    TypeName("cloudBoundaryCollisionNumberFlux");


    // Constructors

        //- Construct from Time and dictionary
        cloudBoundaryCollisionNumberFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudBoundaryCollisionNumberFlux
        (
            const cloudBoundaryCollisionNumberFlux&
        ) = delete;


    //- Destructor
    virtual ~cloudBoundaryCollisionNumberFlux();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudBoundaryCollisionNumberFlux&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
