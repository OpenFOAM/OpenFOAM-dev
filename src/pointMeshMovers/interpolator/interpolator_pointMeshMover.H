/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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
    Foam::pointMeshMovers::interpolator

Description
    Interpolates pointMesh point pointVectorFields

SourceFiles
    interpolator_pointMeshMovers.C

\*---------------------------------------------------------------------------*/

#ifndef interpolator_pointMeshMovers_H
#define interpolator_pointMeshMovers_H

#include "primitiveFields.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class interpolationWeights;

namespace pointMeshMovers
{

/*---------------------------------------------------------------------------*\
                        Class interpolator Declaration
\*---------------------------------------------------------------------------*/

class interpolator
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Name of displacement field
        const word fieldName_;

        const word interpolationScheme_;

        //- Times with pre-specified displacement
        wordList timeNames_;

        //- Times with pre-specified displacement
        scalarField timeVals_;

        //- User-specified interpolator
        mutable autoPtr<interpolationWeights> interpolatorPtr_;

        //- Cached interpolation times
        mutable labelList currentIndices_;

        //- Cached interpolation weights
        mutable scalarField currentWeights_;


public:

    // Constructors

        //- Construct from mesh and dictionary
        interpolator
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~interpolator();


    // Member Functions

        //- Return interpolated pointField for the currentTime
        tmp<pointVectorField> curPointField() const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
