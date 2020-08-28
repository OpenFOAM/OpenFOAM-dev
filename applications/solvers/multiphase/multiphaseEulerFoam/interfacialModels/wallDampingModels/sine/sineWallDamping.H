/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::wallDampingModels::sine

Description

SourceFiles
    sineWallDamping.C

\*---------------------------------------------------------------------------*/

#ifndef sineWallDamping_H
#define sineWallDamping_H

#include "wallDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallDampingModels
{

/*---------------------------------------------------------------------------*\
                            Class sine Declaration
\*---------------------------------------------------------------------------*/

class sine
:
    public wallDampingModel
{
protected:

    // Protected member functions

        //- Return the force limiter field
        virtual tmp<volScalarField> limiter() const;


public:

    //- Runtime type information
    TypeName("sine");


    // Constructors

        //- Construct from components
        sine
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~sine();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallDampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
