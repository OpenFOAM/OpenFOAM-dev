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
    Foam::blendingMethods::stepDisplaced

Description
    Blending method which activates models relating to displacement by a third
    phase based on the value of the phase fraction of that third phase.

    This method requires a volume fractions to be specified for each third
    phase above which non-displaced models are deactivated and displaced models
    are applied instead.

    This model requires a sub-blending method to be specified to determine how
    the dispersed/segregated/general model system (i.e., A_dispersedIn_B,
    B_dispersedIn_A, A_segregatedWith_B, and A_B) is blended. This sub-blending
    is used to generate dispersed/segregated/general coefficients for both the
    non-displaced and displaced regimes, but it will not be used to determine
    the proportion of displaced models that is applied.

    \verbatim
    blending
    {
        default
        {
            type            stepDisplaced;
            displacingAlpha.air none;
            displacingAlpha.water none;
            displacingAlpha.particles 0.15;
            displacingAlpha.bed 0.3;
            blending
            {
                type            linear;
                displaced       no;
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
    }
    \endverbatim

SourceFiles
    stepDisplaced.C

\*---------------------------------------------------------------------------*/

#ifndef stepDisplaced_H
#define stepDisplaced_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                           Class stepDisplaced Declaration
\*---------------------------------------------------------------------------*/

class stepDisplaced
:
    public blendingMethod
{
    // Private Data

        //- The volume fraction of displacing phase above which displaced
        //  models are used
        List<blendingParameter> displacingAlphas_;

        //- List of the displacing phase indices
        List<label> displacingPhases_;

        //- The sub-blending
        autoPtr<blendingMethod> blending_;


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
    TypeName("stepDisplaced");


    // Constructors

        //- Construct from a dictionary and an interface
        stepDisplaced
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~stepDisplaced();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const;

        //- Return whether or not the interface is displaced by a third phase
        virtual bool isDisplacedBy(const label phasei) const;

        //- Return whether this is a function of the volume fractions only
        virtual bool functionOfAlphas() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
