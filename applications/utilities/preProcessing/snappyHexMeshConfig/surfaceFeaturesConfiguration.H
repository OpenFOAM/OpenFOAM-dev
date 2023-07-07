/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::surfaceFeaturesConfiguration

Description
    From a set of input surface geometry files and some configuration
    parameters, writes out a surfacesFeaturesDict configuration file.

SourceFiles
    surfaceFeaturesConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFeaturesConfiguration_H
#define surfaceFeaturesConfiguration_H

#include "caseFileConfiguration.H"
#include "meshingSurfaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class surfaceFeaturesConfiguration Declaration
\*---------------------------------------------------------------------------*/

class surfaceFeaturesConfiguration
:
    public caseFileConfiguration
{
    // Private Data

        //- List of surfaces used for generating features
        const meshingSurfaceList& surfaces_;


    // Private Member Functions

        //- Write surfaces
        void writeSurfaces();

        //- Write subsetFeatures sub-dictionary
        void writeSubsetFeatures();

        //- Write trimFeatures sub-dictionary
        void writeTrimFeatures();

public:

    // Constructors

        //- Construct from components
        surfaceFeaturesConfiguration
        (
            const fileName& name,
            const fileName& dir,
            const Time& time,
            const meshingSurfaceList& surfaces
        );

        //- Disallow default bitwise copy construction
        surfaceFeaturesConfiguration
        (
            const surfaceFeaturesConfiguration&
        ) = delete;


    //- Destructor
    ~surfaceFeaturesConfiguration();


    // Member Functions

        // Writing surfaceFeaturesDict
        void write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceFeaturesConfiguration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
