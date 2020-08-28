/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::solidParticleCloud

Description
    A Cloud of solid particles

SourceFiles
    solidParticleCloudI.H
    solidParticleCloud.C
    solidParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticleCloud_H
#define solidParticleCloud_H

#include "Cloud.H"
#include "solidParticle.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     Class solidParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class solidParticleCloud
:
    public Cloud<solidParticle>
{
    // Private Data

        const fvMesh& mesh_;

        IOdictionary particleProperties_;

        scalar rhop_;
        scalar e_;
        scalar mu_;


public:

    //- Runtime type information
    TypeName("solidParticleCloud");


    // Constructors

        //- Construct given mesh
        solidParticleCloud
        (
            const fvMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Disallow default bitwise copy construction
        solidParticleCloud(const solidParticleCloud&) = delete;


    // Member Functions

        // Access

            inline const fvMesh& mesh() const;

            inline scalar rhop() const;
            inline scalar e() const;
            inline scalar mu() const;


        // Edit

            //- Move the particles under the influence of the given
            //  gravitational acceleration
            void move(const dimensionedVector& g);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidParticleCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
