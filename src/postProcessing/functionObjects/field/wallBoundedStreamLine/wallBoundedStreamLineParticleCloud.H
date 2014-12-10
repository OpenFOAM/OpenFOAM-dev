/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::wallBoundedStreamLineParticleCloud

Description
    A Cloud of streamLine particles

SourceFiles
    streamLineCloud.C

\*---------------------------------------------------------------------------*/

#ifndef streamLineCloud_H
#define streamLineCloud_H

#include "Cloud.H"
#include "wallBoundedStreamLineParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class streamLineCloud Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLineParticleCloud
:
    public Cloud<wallBoundedStreamLineParticle>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallBoundedStreamLineParticleCloud
        (
            const wallBoundedStreamLineParticleCloud&
        );

        //- Disallow default bitwise assignment
        void operator=(const wallBoundedStreamLineParticleCloud&);


public:

    //- Type of parcel the cloud was instantiated for
    typedef wallBoundedStreamLineParticle parcelType;

    // Constructors

        //- Construct given mesh
        wallBoundedStreamLineParticleCloud
        (
            const polyMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Construct from mesh, cloud name, and a list of particles
        wallBoundedStreamLineParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<wallBoundedStreamLineParticle>& particles
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
