/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::reducedUnits

Description

SourceFiles
    reducedUnitsI.H
    reducedUnits.C
    reducedUnitsIO.C

\*---------------------------------------------------------------------------*/

#ifndef reducedUnits_H
#define reducedUnits_H

#include "scalar.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class reducedUnits;

Ostream& operator<<(Ostream&, const reducedUnits&);


/*---------------------------------------------------------------------------*\
                        Class reducedUnits Declaration
\*---------------------------------------------------------------------------*/

class reducedUnits
{
    // Private Data

        // Reduced units

            // Fundamental values

                scalar refLength_;

                scalar refTime_;

                scalar refMass_;

            // Derived values

                scalar refEnergy_;

                scalar refTemp_;

                scalar refForce_;

                scalar refVelocity_;

                scalar refVolume_;

                scalar refPressure_;

                scalar refMassDensity_;

                scalar refNumberDensity_;


    // Private Member Functions

        void calcRefValues();


public:

    // Static Data Members

        //- Static data someStaticData
        static const scalar kb;


    // Constructors

        //- Construct with no argument, uses default values:
        //  length  = 1nm
        //  mass = 1.660538782e-27kg (unified atomic mass unit)
        //  temperature = 1K (therefore, energy = 1*kb)
        reducedUnits();

        //- Construct from components
        reducedUnits
        (
            scalar refLength,
            scalar refTime,
            scalar refMass
        );

        //- Construct from dictionary
        reducedUnits(const IOdictionary& reducedUnitsDict);

        //- Disallow default bitwise copy construction
        reducedUnits(const reducedUnits&) = delete;


    //- Destructor
    ~reducedUnits();


    // Member Functions

        void setRefValues
        (
            scalar refLength,
            scalar refTime,
            scalar refMass
        );

        void setRefValues(const IOdictionary& reducedUnitsDict);


        // Access

            inline scalar refLength() const;

            inline scalar refTime() const;

            inline scalar refMass() const;

            inline scalar refTemp() const;

            inline scalar refEnergy() const;

            inline scalar refForce() const;

            inline scalar refVelocity() const;

            inline scalar refVolume() const;

            inline scalar refPressure() const;

            inline scalar refMassDensity() const;

            inline scalar refNumberDensity() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reducedUnits&) = delete;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const reducedUnits&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reducedUnitsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
