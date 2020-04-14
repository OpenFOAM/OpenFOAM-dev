/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::combustionModels::singleStepCombustion

Description
    Base class for combustion models using multiComponentMixture.

SourceFiles
    singleStepCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef singleStepCombustion_H
#define singleStepCombustion_H

#include "multiComponentMixture.H"
#include "ThermoCombustion.H"
#include "ReactionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                    Class singleStepCombustion Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class singleStepCombustion
:
    public ThermoCombustion<ReactionThermo>
{
protected:

    // Protected data

        const multiComponentMixture<ThermoType>& mixture_;

        //- The single-step reaction
        autoPtr<Reaction<ThermoType>> reaction_;

        //- Stoichiometric air-fuel mass ratio
        dimensionedScalar stoicRatio_;

        //- Stoichiometric oxygen-fuel mass ratio
        dimensionedScalar s_;

        //- Heat of combustion [J/Kg]
        dimensionedScalar qFuel_;

        //- Stoichiometric coefficient for the reaction.
        scalarList specieStoichCoeffs_;

        //- Mass concentrations at stoichiometric mixture for fres.
        scalarList Yprod0_;

        //- List of components residual
        PtrList<volScalarField> fres_;

        //- Inert specie index
        label inertIndex_;

        //- Fuel specie index
        label fuelIndex_;

        //- List to indicate if specie is produced/consumed
        List<int> specieProd_;

        //- Fuel consumption rate
        volScalarField wFuel_;

        //- Semi-implicit (true) or explicit (false) treatment
        bool semiImplicit_;


    // Protected member functions

        //- Calculate qFuel
        void calculateqFuel();

        //- Calculate air/fuel and oxygen/fuel ratio
        void massAndAirStoichRatios();

        //- Calculate maximum products at stoichiometric mixture
        void calculateMaxProducts();


public:

    // Constructors

        //- Construct from components
        singleStepCombustion
        (
            const word& modelType,
            const ReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        singleStepCombustion(const singleStepCombustion&);


    //- Destructor
    virtual ~singleStepCombustion();


    // Member Functions

        // Access functions

    const multiComponentMixture<ThermoType>& mixture() const
    {
        return mixture_;
    }

            //- Return the single step reaction
            inline const Reaction<ThermoType>& reaction() const;

            //- Return the stoichiometric air-fuel mass ratio
            inline const dimensionedScalar stoicRatio() const;

            //- Return the Stoichiometric oxygen-fuel mass ratio
            inline const dimensionedScalar s() const;

            //- Return the heat of combustion [J/Kg]
            inline const dimensionedScalar qFuel() const;

            //- Return the stoichiometric coefficient for the reaction
            inline const List<scalar>& specieStoichCoeffs() const;

            //- Return the list of components residual
            inline tmp<volScalarField> fres(const label index) const;

            //- Return the inert specie index
            inline label inertIndex() const;

            //- Return the fuel specie index
            inline label fuelIndex() const;

            //- Return the list to indicate if specie is produced/consumed
            inline const List<int>& specieProd() const;

            //- Return the list of products mass concentrations
            inline const scalarList& Yprod0() const;


        //- Fuel consumption rate matrix
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();

        //- Calculates the residual for all components
        void fresCorrect();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleStepCombustion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "singleStepCombustionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleStepCombustion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
