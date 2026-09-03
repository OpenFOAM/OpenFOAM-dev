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
    Foam::blendingMethods::zonal

Description
    Blending method which applies different sub-blending methods to different
    cell zones. The blending methods are specified in a list of
    sub-dictionaries, with the cell zone defined by the sub-dictionary names. A
    default method is also specified for any part of the mesh not within one of
    the identified zones.

    Where multiple zones overlap then coefficients relating to latter zones
    will take precedence.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            zonal;
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
            zones
            {
                myCellZone0
                {
                    type            segregated;
                }
                myCellZone1
                {
                    type            none;
                }
            }
        }
    }
    \endverbatim

SourceFiles
    zonal.C
    zonalBoundaryIndices.C

\*---------------------------------------------------------------------------*/

#ifndef zonal_H
#define zonal_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                           Class zonal Declaration
\*---------------------------------------------------------------------------*/

class zonal
:
    public blendingMethod
{
    // Private classes

        //- Zone boundary indexing information
        class boundaryIndices;


    // Private Data

        //- The blending method outside the zones
        autoPtr<blendingMethod> defaultBlending_;

        //- The indices of the zones for which non-default blending applies
        labelList zoneIndices_;

        //- The blending methods in the zones
        PtrList<blendingMethod> zoneBlendings_;


    // Private Member Functions

        //- Combine default and zone-specific function values
        tmp<volScalarField> f
        (
            const tmp<volScalarField>& defaultF,
            const PtrList<volScalarField>& zoneFs
        ) const;


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
    TypeName("zonal");


    // Constructors

        //- Construct from a dictionary and an interface
        zonal
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~zonal();


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
