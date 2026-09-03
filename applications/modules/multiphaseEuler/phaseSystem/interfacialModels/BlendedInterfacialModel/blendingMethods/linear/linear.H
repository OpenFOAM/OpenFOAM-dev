/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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
    considered to become continuous to be specified for each phase.

    Alternatively, these volume fractions can be omitted or replaced with the
    keyword 'none' to represent a phase which cannot become continuous. E.g.,
    a particulate phase.

    If there are more than two phases, this model requires an additional
    'displaced' switch. If 'displaced' is set to false, then displaced models
    will not be used. If it is set to 'true', then displaced models will be
    activated as the combined phase fraction of the two phases to which the
    model applies reduces relative to the other phases. Displaced and
    non-displaced models will be combined using the same linear function as
    is used to blend the dispersed and segregated models.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            linear;
            displaced       yes;
            minFullyContinuousAlpha.air 0.8;
            minPartlyContinuousAlpha.air 0.5;
            minFullyContinuousAlpha.water 0.8;
            minPartlyContinuousAlpha.water 0.5;
            minFullyContinuousAlpha.particles none;
            minPartlyContinuousAlpha.particles none;
            minFullyContinuousAlpha.bed none;
            minPartlyContinuousAlpha.bed none;
        }
    }
    \endverbatim

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

        //- Apply displaced models? The displaced model coefficient will be
        //  calculated as one-minus the continuous coefficient calculated for
        //  the two-phase sub-system.
        const bool displaced_;

        //- Minimum fraction of phases which can be considered fully continuous
        Pair<blendingParameter> minFullyContinuousAlpha_;

        //- Minimum fraction of phases which can be considered partly continuous
        Pair<blendingParameter> minPartlyContinuousAlpha_;


protected:

    // Protected Member Functions

        //- Return the coefficient for models for which the indexed phase is
        //  considered continuous
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label index
        ) const;

        //- Return the coefficient for when the interface is displaced by a
        //  third phase
        virtual tmp<volScalarField> fDisplaced
        (
            const UPtrList<const volScalarField>& alphas,
            const label displacingPhasei
        ) const;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from a dictionary and an interface
        linear
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~linear();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const;

        //- Return whether or not the interface is displaced by a third phase
        virtual bool isDisplacedBy(const label displacingPhasei) const;

        //- Return whether this is a function of the volume fractions only
        virtual bool functionOfAlphas() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
