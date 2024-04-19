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
    Foam::LandauTellerReactionRate

Description
    Landau-Teller reaction rate.

SourceFiles
    LandauTellerReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef LandauTellerReactionRate_H
#define LandauTellerReactionRate_H

#include "speciesTable.H"
#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class LandauTellerReactionRate;

Ostream& operator<<(Ostream&, const LandauTellerReactionRate&);


/*---------------------------------------------------------------------------*\
                  Class LandauTellerReactionRate Declaration
\*---------------------------------------------------------------------------*/

class LandauTellerReactionRate
{
    // Private Data

        scalar beta_;
        scalar A_;
        scalar Ta_;
        scalar B_;
        scalar C_;


public:

    // Constructors

        //- Construct from components
        inline LandauTellerReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const scalar B,
            const scalar C
        );

        //- Construct from dictionary
        inline LandauTellerReactionRate
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "LandauTeller";
        }

        //- Pre-evaluation hook
        inline void preEvaluate() const;

        //- Post-evaluation hook
        inline void postEvaluate() const;

        //- Return the rate
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- The derivative of the rate w.r.t. temperature
        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Is the rate a function of concentration?
        inline bool hasDdc() const;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& ddc
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const LandauTellerReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LandauTellerReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
