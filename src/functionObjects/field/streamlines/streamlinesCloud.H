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
    Foam::streamlinesCloud

Description
    A Cloud of streamlines particles

SourceFiles
    streamlinesCloud.C

\*---------------------------------------------------------------------------*/

#ifndef streamlinesCloud_H
#define streamlinesCloud_H

#include "Cloud.H"
#include "streamlinesParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class streamlinesCloud Declaration
\*---------------------------------------------------------------------------*/

class streamlinesCloud
:
    public Cloud<streamlinesParticle>
{
public:

    //- Type of parcel the cloud was instantiated for
    typedef streamlinesParticle parcelType;


    //- Runtime type information
    TypeName("streamlinesCloud");


    // Constructors

        //- Construct given mesh
        streamlinesCloud
        (
            const polyMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Construct from mesh, cloud name, and a list of particles
        streamlinesCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<streamlinesParticle>& particles
        );

        //- Disallow default bitwise copy construction
        streamlinesCloud(const streamlinesCloud&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const streamlinesCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
