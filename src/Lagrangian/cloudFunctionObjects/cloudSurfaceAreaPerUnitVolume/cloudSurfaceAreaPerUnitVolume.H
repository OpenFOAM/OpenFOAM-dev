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
    Foam::functionObjects::cloudSurfaceAreaPerUnitVolume

Description
    Function to write the surface area per unit volume field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudSurfaceAreaPerUnitVolume1
    {
        type            cloudSurfaceAreaPerUnitVolume;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudSurfaceAreaPerUnitVolume.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudSurfaceAreaPerUnitVolume_H
#define functionObjects_cloudSurfaceAreaPerUnitVolume_H

#include "cloudFvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
              Class cloudSurfaceAreaPerUnitVolume Declaration
\*---------------------------------------------------------------------------*/

class cloudSurfaceAreaPerUnitVolume
:
    public cloudFvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("cloudSurfaceAreaPerUnitVolume");


    // Public Static Data

        //- The name of the surface area per unit volume field
        const word AvName_;


    // Constructors

        //- Construct from Time and dictionary
        cloudSurfaceAreaPerUnitVolume
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudSurfaceAreaPerUnitVolume
        (
            const cloudSurfaceAreaPerUnitVolume&
        ) = delete;


    //- Destructor
    virtual ~cloudSurfaceAreaPerUnitVolume();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate and store the surface area per unit volume
        virtual bool execute();

        //- Write the surface area per unit volume
        virtual bool write();

        //- Clear the surface area per unit volume
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudSurfaceAreaPerUnitVolume&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
