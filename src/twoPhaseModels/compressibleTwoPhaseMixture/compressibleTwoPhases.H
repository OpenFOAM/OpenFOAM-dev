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
    Foam::compressibleTwoPhaseMixture

Description
    Interface to two rhoThermo-based phases

SourceFiles
    compressibleTwoPhases.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTwoPhases_H
#define compressibleTwoPhases_H

#include "twoPhases.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class compressibleTwoPhases Declaration
\*---------------------------------------------------------------------------*/

class compressibleTwoPhases
:
    virtual public twoPhases
{
public:

    //- Runtime type information
    TypeName("compressibleTwoPhases");


    // Constructors

        //- Default construct
        compressibleTwoPhases()
        {}


    //- Destructor
    virtual ~compressibleTwoPhases()
    {}


    // Member Functions

        //- Return the thermo for phase 1
        virtual const rhoThermo& thermo1() const = 0;

        //- Return the thermo for phase 2
        virtual const rhoThermo& thermo2() const = 0;

        //- Return the density of a given phase
        inline const rhoThermo& thermo(const bool index) const
        {
            return index ? thermo2() : thermo1();
        }

        //- Return the density of phase 1
        virtual const volScalarField& rho1() const = 0;

        //- Return the density of phase 2
        virtual const volScalarField& rho2() const = 0;

        //- Return the density of a given phase
        inline const volScalarField& rho(const bool index) const
        {
            return index ? rho2() : rho1();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
