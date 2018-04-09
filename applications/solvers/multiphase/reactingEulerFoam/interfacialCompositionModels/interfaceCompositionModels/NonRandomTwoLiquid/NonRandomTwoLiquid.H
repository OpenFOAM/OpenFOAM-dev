/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::NonRandomTwoLiquid

Description
    Non ideal law for the mixing of two species. A separate composition model
    is given for each species. The composition of a species is equal to the
    value given by the model, scaled by the species fraction in the bulk of the
    other phase, and multiplied by the activity coefficient for that species.
    The gas behaviour is assumed ideal; i.e. the fugacity coefficient is taken
    as equal to 1.

SourceFiles
    NonRandomTwoLiquid.C

\*---------------------------------------------------------------------------*/

#ifndef NonRandomTwoLiquid_H
#define NonRandomTwoLiquid_H

#include "InterfaceCompositionModel.H"
#include "saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace interfaceCompositionModels
{

/*---------------------------------------------------------------------------*\
                         Class NonRandomTwoLiquid Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class NonRandomTwoLiquid
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Activity coefficient for species 1
        volScalarField gamma1_;

        //- Activity coefficient for species 2
        volScalarField gamma2_;

        //- Name of species 1
        word species1Name_;

        //- Name of species 2
        word species2Name_;

        //- Indiex of species 1 within this thermo
        label species1Index_;

        //- Indiex of species 2 within this thermo
        label species2Index_;

        //- Non-randomness constant parameter for species 1
        dimensionedScalar alpha12_;

        //- Non-randomness constant parameter for species 2
        dimensionedScalar alpha21_;

        //- Non-randomness linear parameter for species 1
        dimensionedScalar beta12_;

        //- Non-randomness linear parameter for species 2
        dimensionedScalar beta21_;

        //- Interaction parameter model for species 1
        autoPtr<saturationModel> saturationModel12_;

        //- Interaction parameter model for species 2
        autoPtr<saturationModel> saturationModel21_;

        //- Composition model for species 1
        autoPtr<interfaceCompositionModel> speciesModel1_;

        //- Composition model for species 2
        autoPtr<interfaceCompositionModel> speciesModel2_;


public:

    //- Runtime type information
    TypeName("nonRandomTwoLiquid");

    // Constructors

        //- Construct from components
        NonRandomTwoLiquid
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~NonRandomTwoLiquid();


    // Member Functions

        //- Update the composition
        virtual void update(const volScalarField& Tf);

        //- The interface species fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- The interface species fraction derivative w.r.t. temperature
        virtual tmp<volScalarField> YfPrime
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfaceCompositionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NonRandomTwoLiquid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
