/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::LESModels::cubeRootVolDelta

Description
    Simple cube-root of cell volume delta used in LES SGS models.

SourceFiles
    cubeRootVolDelta.C

\*---------------------------------------------------------------------------*/

#ifndef cubeRootVolDelta_H
#define cubeRootVolDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class cubeRootVolDelta Declaration
\*---------------------------------------------------------------------------*/

class cubeRootVolDelta
:
    public LESdelta
{
    // Private Data

        scalar deltaCoeff_;


    // Private Member Functions

        // Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("cubeRootVol");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        cubeRootVolDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        cubeRootVolDelta(const cubeRootVolDelta&) = delete;


    //- Destructor
    virtual ~cubeRootVolDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cubeRootVolDelta&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
