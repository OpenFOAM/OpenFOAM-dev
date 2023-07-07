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
    Foam::blockMeshConfigurationBase

Description
    Functions to configure and write a blockMeshDict configuration file.

SourceFiles
    blockMeshConfigurationBase.C

\*---------------------------------------------------------------------------*/

#ifndef blockMeshConfigurationBase_H
#define blockMeshConfigurationBase_H

#include "caseFileConfiguration.H"
#include "meshingSurfaceList.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class blockMeshConfigurationBase Declaration
\*---------------------------------------------------------------------------*/

class blockMeshConfigurationBase
:
    public caseFileConfiguration
{
protected:

    // Protected Data

        //- Bounding box for the background mesh block
        boundBox bb_;

        //- Optional settings for patch names and types in the background mesh
        HashTable<Pair<word>> patchOpts_;


    // Protected Member Functions

        //- Round a bounding box by the rounding scale
        void roundBoundingBox(boundBox& bb, const scalar s);

        //- Parse the patch commandline options
        Pair<word> readPatchOption(const word& option) const;

        //- Write a vertex entry in the vertices list
        void writeVertex
        (
            const word& x,
            const word& y,
            const word& z
        );


public:

    // Constructors

        //- Construct from components
        blockMeshConfigurationBase
        (
            const fileName& name,
            const fileName& dir,
            const Time& time,
            const meshingSurfaceList& surfaces,
            const HashTable<Pair<word>>& patchOpts
        );

        //- Disallow default bitwise copy construction
        blockMeshConfigurationBase(const blockMeshConfigurationBase&) = delete;


    //- Destructor
    ~blockMeshConfigurationBase();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blockMeshConfigurationBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
