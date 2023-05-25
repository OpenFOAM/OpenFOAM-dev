/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::blendingMethods::linear

Description
    Blending method based on piecewise linear functions. Supports the full
    range of phase fraction space. E.g., from droplets in air, through a
    segregated regime, to bubbly flow.

    This method requires two volume fractions between which the phase is
    considered to become continuous to be specified for both phases.

    Alternatively, these volume fractions can be omitted or replaced with the
    keyword "none" to represent a phase which cannot become continuous. E.g.,
    a particulate phase.

SourceFiles
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

class linear
:
    public blendingMethod
{
    // Private Data

        //- Minimum fraction of phases which can be considered fully continuous
        Pair<blendingParameter> minFullyContinuousAlpha_;

        //- Minimum fraction of phases which can be considered partly continuous
        Pair<blendingParameter> minPartlyContinuousAlpha_;


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
    TypeName("linear");


    // Constructors

        //- Construct from a dictionary and an interface
        linear
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    ~linear();


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
