/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::incompressibleTwoPhases

Description
    Interface to two constant density phases

SourceFiles
    incompressibleTwoPhases.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhases_H
#define incompressibleTwoPhases_H

#include "twoPhases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class incompressibleTwoPhases Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhases
:
    virtual public twoPhases
{
public:

    //- Runtime type information
    TypeName("incompressibleTwoPhases");


    // Constructors

        //- Default construct
        incompressibleTwoPhases()
        {}


    //- Destructor
    virtual ~incompressibleTwoPhases()
    {}


    // Member Functions

        //- Return the density of phase1
        virtual const dimensionedScalar& rho1() const = 0;

        //- Return the density of phase2
        virtual const dimensionedScalar& rho2() const = 0;

        //- Return the density of a given phase
        inline const dimensionedScalar& rho(const bool index) const
        {
            return index ? rho2() : rho1();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
