/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::immiscibleCompressibleTwoPhaseMixture

Description
    Class to represent a mixture of two immiscible rhoThermo-based phases

SourceFiles
    immiscibleCompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleCompressibleTwoPhaseMixture_H
#define immiscibleCompressibleTwoPhaseMixture_H

#include "compressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class immiscibleCompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleCompressibleTwoPhaseMixture
:
    public compressibleTwoPhaseMixture,
    public interfaceProperties
{

public:

    // Constructors

        //- Construct from the mixture velocity
        immiscibleCompressibleTwoPhaseMixture(const volVectorField& U);


    //- Destructor
    virtual ~immiscibleCompressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Correct the mixture and interface properties
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
