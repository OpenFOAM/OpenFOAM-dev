/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::blendingMethods::hyperbolic

Description

SourceFiles
    hyperbolic.C

\*---------------------------------------------------------------------------*/

#ifndef hyperbolic_H
#define hyperbolic_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class hyperbolic Declaration
\*---------------------------------------------------------------------------*/

class hyperbolic
:
    public blendingMethod
{
    // Private Data

        //- Maximum fraction of phases which can be considered dispersed
        HashTable<dimensionedScalar, word, word::hash> maxDispersedAlpha_;

        //- Width of the transition
        const dimensionedScalar transitionAlphaScale_;


public:

    //- Runtime type information
    TypeName("hyperbolic");


    // Constructors

        //- Construct from a dictionary and a list of phase names
        hyperbolic
        (
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    ~hyperbolic();


    // Member Functions

        //- Factor for first phase
        virtual tmp<volScalarField> f1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for second phase
        virtual tmp<volScalarField> f2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
