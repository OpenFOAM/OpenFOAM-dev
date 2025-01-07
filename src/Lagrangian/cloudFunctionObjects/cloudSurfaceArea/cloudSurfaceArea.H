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
    Foam::functionObjects::cloudSurfaceArea

Description
    Function to write the surface area field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudSurfaceArea1
    {
        type            cloudSurfaceArea;
        libs            ("libfieldLagrangianFunctionObjects.so");
        Lagrangian      cloud;
        field           m;
        weightField     number;
    }
    \endverbatim

SourceFiles
    cloudSurfaceArea.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudSurfaceArea_H
#define functionObjects_cloudSurfaceArea_H

#include "cloudLagrangianMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudSurfaceArea Declaration
\*---------------------------------------------------------------------------*/

class cloudSurfaceArea
:
    public cloudLagrangianMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("cloudSurfaceArea");


    // Public Static Data

        //- The name of the surface area field
        static const word aName_;


    // Constructors

        //- Construct from Time and dictionary
        cloudSurfaceArea
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudSurfaceArea(const cloudSurfaceArea&) = delete;


    //- Destructor
    virtual ~cloudSurfaceArea();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate and store the surface areas
        virtual bool execute();

        //- Update the surface areas
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Write the surface areas
        virtual bool write();

        //- Clear the surface areas
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudSurfaceArea&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
