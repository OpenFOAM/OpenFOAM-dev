/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::sampledSetCloud

Description
    A Cloud of sampledSet particles

SourceFiles
    sampledSetCloud.C

\*---------------------------------------------------------------------------*/

#ifndef sampledSetCloud_H
#define sampledSetCloud_H

#include "Cloud.H"
#include "sampledSetParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sampledSetCloud Declaration
\*---------------------------------------------------------------------------*/

class sampledSetCloud
:
    public Cloud<sampledSetParticle>
{
public:

    //- Type of parcel the cloud was instantiated for
    typedef sampledSetParticle parcelType;


    //- Runtime type information
    TypeName("sampledSetCloud");


    // Constructors

        //- Construct given mesh
        sampledSetCloud
        (
            const polyMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Construct from mesh, cloud name, and a list of particles
        sampledSetCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<sampledSetParticle>& particles
        );

        //- Disallow default bitwise copy construction
        sampledSetCloud(const sampledSetCloud&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sampledSetCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
