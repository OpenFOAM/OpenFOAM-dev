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
    Foam::blendingMethods::displaced

Description
    Blending method for the case in which the models relating to displacement
    by a third phase should always be used. This is unlikely to be used on its
    own. It is intended that this be used within another blending method which
    applies it to the subset of the mesh in which the third phase is present.

    This model requires a sub-blending method to be specified to determine how
    the different displaced models (i.e., A_dispersedIn_B_displacedBy_C,
    B_dispersedIn_A_displacedBy_C, A_segregatedWith_B_displacedBy_C, and
    A_B_displacedBy_C) are to be combined. 'none' is a valid and likely choice
    for this sub-blending model, and will result in only the default/general
    model (i.e., A_B_displacedBy_C) being used.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            zonal;

            // Use linear blending outside the 'bed' zone
            default
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

            // Use displaced blending inside the bed zone, and
            // only apply the default/general displaced model
            zones
            {
                bed
                {
                    type            displaced;
                    phase           bed;
                    blending
                    {
                        type            none;
                    }
                }
            }
        }
    }
    \endverbatim

SourceFiles
    displaced.C

\*---------------------------------------------------------------------------*/

#ifndef displaced_H
#define displaced_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                           Class displaced Declaration
\*---------------------------------------------------------------------------*/

class displaced
:
    public blendingMethod
{
protected:

    // Private Data

        //- The index of the displacing phase
        const label phasei_;

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
    TypeName("displaced");


    // Constructors

        //- Construct from a dictionary and an interface
        displaced
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~displaced();


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
