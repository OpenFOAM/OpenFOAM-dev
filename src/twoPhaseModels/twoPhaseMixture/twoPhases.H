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
    Foam::twoPhases

Description
    Interface to two phases

SourceFiles
    twoPhases.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhases_H
#define twoPhases_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class twoPhases Declaration
\*---------------------------------------------------------------------------*/

class twoPhases
{
public:

    TypeName("twoPhases");


    // Constructors

        //- Default construct
        twoPhases()
        {}


    //- Destructor
    virtual ~twoPhases()
    {}


    // Member Functions

        //- Access the mesh
        const fvMesh& mesh() const
        {
            return alpha1().mesh();
        }

        //- Return the phase-fraction of phase 1
        virtual const volScalarField& alpha1() const = 0;

        //- Return the phase-fraction of phase 2
        virtual const volScalarField& alpha2() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
