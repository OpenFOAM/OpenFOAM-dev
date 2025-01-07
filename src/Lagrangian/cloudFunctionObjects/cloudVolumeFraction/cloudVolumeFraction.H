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
    Foam::functionObjects::cloudVolumeFraction

Description
    Function to write the volume fraction field for a cloud

Usage
    \table
        Property      | Description                  | Required? | Default
        cloud         | Name of the cloud            | yes       |
    \endtable

    Example specification:
    \verbatim
    cloudVolumeFraction1
    {
        type            cloudVolumeFraction;
        libs            ("libLagrangianCloudFunctionObjects.so");
        cloud           cloud;
    }
    \endverbatim

SourceFiles
    cloudVolumeFraction.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudVolumeFraction_H
#define functionObjects_cloudVolumeFraction_H

#include "cloudFvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudVolumeFraction Declaration
\*---------------------------------------------------------------------------*/

class cloudVolumeFraction
:
    public cloudFvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("cloudVolumeFraction");


    // Public Static Data

        //- The name of the volume fraction field
        const word alphaName_;


    // Constructors

        //- Construct from Time and dictionary
        cloudVolumeFraction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudVolumeFraction(const cloudVolumeFraction&) = delete;


    //- Destructor
    virtual ~cloudVolumeFraction();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate and store the volume fraction
        virtual bool execute();

        //- Write the volume fraction
        virtual bool write();

        //- Clear the volume fraction
        virtual bool clear();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudVolumeFraction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
