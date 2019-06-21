/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::singleStepReactingMixture

Description
    Single step reacting mixture

SourceFiles
    singleStepReactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef singleStepReactingMixture_H
#define singleStepReactingMixture_H

#include "chemistryReader.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class singleStepReactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class singleStepReactingMixture
:
    public reactingMixture<ThermoType>
{

protected:

    // Protected data

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


    // Protected member functions

        //- Calculate qFuel
        void calculateqFuel();

        //- Calculate air/fuel and oxygen/fuel ratio
        void massAndAirStoichRatios();

        //- Calculate maximum products at stoichiometric mixture
        void calculateMaxProducts();


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        singleStepReactingMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        singleStepReactingMixture(const singleStepReactingMixture&) = delete;


    //- Destructor
    virtual ~singleStepReactingMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "singleStepReactingMixture<" + ThermoType::typeName() + '>';
        }

        //- Calculates the residual for all components
        void fresCorrect();


        // Access functions

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


        // I-O

            //- Read dictionary
            void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singleStepReactingMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "singleStepReactingMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "singleStepReactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
