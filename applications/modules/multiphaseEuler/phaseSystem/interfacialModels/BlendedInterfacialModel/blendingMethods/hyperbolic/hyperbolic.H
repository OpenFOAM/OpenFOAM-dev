/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2025 OpenFOAM Foundation
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
    Blending method based on smooth hyperbolic functions. Supports the full
    range of phase fraction space. E.g., from droplets in air, through a
    segregated regime, to bubbly flow.

    This method requires a volume fraction at which the phase is considered to
    become continuous, and a range of volume fraction over which this occurs,
    to be specified for both phases.

    Alternatively, the volume fraction can be omitted or replaced with the
    keyword "none" to represent a phase which cannot become continuous. E.g.,
    a particulate phase.

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

        //- Minimum fraction of phases which can be considered continuous
        Pair<blendingParameter> minContinuousAlpha_;

        //- Width of the transition
        const blendingParameter transitionAlphaScale_;


protected:

    // Protected Member Functions

        //- Evaluate the blending function
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label phaseSet,
            const label systemSet
        ) const;


public:

    //- Runtime type information
    TypeName("hyperbolic");


    // Constructors

        //- Construct from a dictionary and an interface
        hyperbolic
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    ~hyperbolic();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
