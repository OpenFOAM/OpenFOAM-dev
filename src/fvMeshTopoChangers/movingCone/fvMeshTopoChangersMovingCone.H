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
    Foam::fvMeshTopoChangers::movingCone

Description
    Sample fvMeshTopoChanger that moves an object in x direction and
    introduces/removes layers.

SourceFiles
    fvMeshTopoChangersMovingCone.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersMovingCone_H
#define fvMeshTopoChangersMovingCone_H

#include "fvMeshTopoChanger.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                    Class movingCone Declaration
\*---------------------------------------------------------------------------*/

class movingCone
:
    public fvMeshTopoChanger
{
    // Private Data

        polyTopoChanger topoChanger_;

        //- Motion dictionary
        dictionary motionDict_;

        //- Motion velocity amplitude
        vector motionVelAmplitude_;

        //- Motion velocity period
        scalar motionVelPeriod_;

        //- Motion velocity period
        vector curMotionVel_;

        //- Left edge
        scalar leftEdge_;

        //- Current left obstacle position
        scalar curLeft_;

        //- Current right obstacle position
        scalar curRight_;

        //- Vertex motion mask
        scalarField motionMask_;


    // Private Member Functions

        //- Add mixer zones and modifiers
        void addZonesAndModifiers();

        //- Markup motion vertices
        tmp<scalarField> vertexMarkup
        (
            const pointField& p,
            const scalar curLeft,
            const scalar curRight
        ) const;


public:

    //- Runtime type information
    TypeName("movingCone");


    // Constructors

        //- Construct from fvMesh
        explicit movingCone(fvMesh& io);

        //- Disallow default bitwise copy construction
        movingCone(const movingCone&) = delete;


    //- Destructor
    virtual ~movingCone();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void updateMesh(const mapPolyMesh&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const mapDistributePolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const movingCone&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
