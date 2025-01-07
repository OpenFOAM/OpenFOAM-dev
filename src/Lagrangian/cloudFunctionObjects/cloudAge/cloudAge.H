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
    Foam::functionObjects::cloudAge

Description
    Function to compute and write an age field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudAge1
    {
        type            cloudAge;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudAge.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudAge_H
#define functionObjects_cloudAge_H

#include "cloudLagrangianMeshFunctionObject.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudAge Declaration
\*---------------------------------------------------------------------------*/

class cloudAge
:
    public cloudLagrangianMeshFunctionObject
{
private:

    // Private Data

        //- Age field
        LagrangianScalarDynamicField age_;


public:

    //- Runtime type information
    TypeName("cloudAge");


    // Constructors

        //- Construct from Time and dictionary
        cloudAge
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudAge(const cloudAge&) = delete;


    //- Destructor
    virtual ~cloudAge();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate and store the volumes
        virtual bool execute();

        //- Update the ages
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Write the volumes
        virtual bool write();

        //- Clear the volumes
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudAge&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
