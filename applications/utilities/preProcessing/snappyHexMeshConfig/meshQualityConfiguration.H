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
    Foam::meshQualityConfiguration

Description
    Writes a meshQualityDict file which is included from the snappyHexMeshDict
    file.

SourceFiles
    meshQualityConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef meshQualityConfiguration_H
#define meshQualityConfiguration_H

#include "caseFileConfiguration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class meshQualityConfiguration Declaration
\*---------------------------------------------------------------------------*/

class meshQualityConfiguration
:
    public caseFileConfiguration
{

public:

    // Constructors

        meshQualityConfiguration
        (
            const fileName& name,
            const fileName& dir,
            const Time& time
        );

        //- Disallow default bitwise copy construction
        meshQualityConfiguration(const meshQualityConfiguration&) = delete;

    //- Destructor
    ~meshQualityConfiguration();


    // Member Functions

        // Write the snappyHexMeshDict
        void write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshQualityConfiguration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
