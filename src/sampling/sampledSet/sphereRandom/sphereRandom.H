/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::sampledSets::sphereRandom

Description
    Random samples within a sphere

Usage
    \table
        Property    | Description                            | Req'd? | Default
        centre      | Centre of the sphere                   | yes    |
        radius      | Radius of the sphere                   | yes    |
        nPoints     | The number of points                   | yes    |
        axis        | The coordinate axis that is written    | yes    |
    \endtable

    Example specification:
    \verbatim
    {
        type        sphereRandom;
        centre      (0.95 0 0.25);
        radius      0.25;
        nPoints     200;
        axis        x;
    }
    \endverbatim

SourceFiles
    sphereRandom.C

\*---------------------------------------------------------------------------*/

#ifndef sphereRandom_H
#define sphereRandom_H

#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSets
{

/*---------------------------------------------------------------------------*\
                        Class sphereRandom Declaration
\*---------------------------------------------------------------------------*/

class sphereRandom
:
    public sampledSet
{
    // Private Data

        //- Centre point
        const point centre_;

        //- Radius
        const scalar radius_;

        //- Number of points
        const label nPoints_;


    // Private Member Functions

        //- Calculate all the sampling points
        virtual bool calcSamples
        (
            DynamicList<point>& samplingPositions,
            DynamicList<scalar>&,
            DynamicList<label>& samplingSegments,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces
        ) const;


public:

    //- Runtime type information
    TypeName("sphereRandom");


    // Constructors

        //- Construct from dictionary
        sphereRandom
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sphereRandom();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sampledSets
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
