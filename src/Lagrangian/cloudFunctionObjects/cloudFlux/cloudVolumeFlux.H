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
    Foam::functionObjects::cloudVolumeFlux

Description
    Function to write the number flux field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudVolumeFlux1
    {
        type            cloudVolumeFlux;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudVolumeFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudVolumeFlux_H
#define functionObjects_cloudVolumeFlux_H

#include "cloudFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudVolumeFlux Declaration
\*---------------------------------------------------------------------------*/

class cloudVolumeFlux
:
    public cloudFlux
{
private:

    // Private Member Functions

        //- Return the property for which to compute the flux
        virtual tmp<LagrangianSubScalarSubField> q
        (
            const LagrangianSubScalarSubField& fraction
        ) const;


public:

    //- Runtime type information
    TypeName("cloudVolumeFlux");


    // Constructors

        //- Construct from Time and dictionary
        cloudVolumeFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudVolumeFlux(const cloudVolumeFlux&) = delete;


    //- Destructor
    virtual ~cloudVolumeFlux();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudVolumeFlux&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
