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
    Foam::functionObjects::cloudVolume

Description
    Function to write the volume field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudVolume1
    {
        type            cloudVolume;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudVolume.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudVolume_H
#define functionObjects_cloudVolume_H

#include "cloudLagrangianMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudVolume Declaration
\*---------------------------------------------------------------------------*/

class cloudVolume
:
    public cloudLagrangianMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("cloudVolume");


    // Constructors

        //- Construct from Time and dictionary
        cloudVolume
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudVolume(const cloudVolume&) = delete;


    //- Destructor
    virtual ~cloudVolume();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate and store the volumes
        virtual bool execute();

        //- Update the volumes
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
        void operator=(const cloudVolume&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
