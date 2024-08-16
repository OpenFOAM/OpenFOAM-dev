/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::inhomogeneousEGRMixture

Description
    Combustion mixture for inhomogeneous lean, stoichiometric or rich conditions
    with exhaust gas recirculation

    in which the mixture-fraction \c ft, fuel mass-fraction \c fu and exhaust
    gas recirculation mass-fraction \c egr define the composition.

SourceFiles
    inhomogeneousEGRMixture.C

\*---------------------------------------------------------------------------*/

#ifndef inhomogeneousEGRMixture_H
#define inhomogeneousEGRMixture_H

#include "dimensionedTypes.H"
#include "FieldListSlice.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class inhomogeneousEGRMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class inhomogeneousEGRMixture
{
public:

    // Public Typedefs

        //- The type of thermodynamics this mixture is instantiated for
        typedef ThermoType thermoType;

        //- Mixing type for thermodynamic properties
        typedef ThermoType thermoMixtureType;

        //- Mixing type for transport properties
        typedef ThermoType transportMixtureType;


private:

    // Private Enumerations

        //- Specie indices
        enum species { FT, FU, B, EGR };


    // Private Data

        //- Stoichiometric air/fuel ratio
        dimensionedScalar stoicRatio_;

        //- Fuel thermodynamic model
        thermoType fuel_;

        //- Oxidant thermodynamic model
        thermoType oxidant_;

        //- Product thermodynamic model
        thermoType products_;

        //- Mutable storage for the mixed thermodynamic model
        mutable thermoType mixture_;


public:

    // Constructors

        //- Construct from a dictionary
        inhomogeneousEGRMixture(const dictionary&);

        //- Disallow default bitwise copy construction
        inhomogeneousEGRMixture
        (
            const inhomogeneousEGRMixture<ThermoType>&
        ) = delete;


    //- Destructor
    virtual ~inhomogeneousEGRMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "inhomogeneousEGRMixture<" + ThermoType::typeName() + '>';
        }

        //- Return the specie names
        static wordList specieNames()
        {
            return {"ft", "fu", "b", "egr"};
        }

        //- Return the residual fraction of fuel in the burnt mixture
        scalar fres(const scalarFieldListSlice& Y) const;

        //- Return the mixture for the given composition
        const thermoType& mixture
        (
            const scalar ft,
            const scalar fu,
            const scalar egr
        ) const;

        //- Return the mixture for thermodynamic properties
        const thermoMixtureType& thermoMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&
        ) const;

        //- Return the mixture for transport properties
        const transportMixtureType& transportMixture
        (
            const scalarFieldListSlice&,
            const thermoMixtureType&
        ) const;

        //- Return the reactant mixture
        const thermoType& reactants(const scalarFieldListSlice&) const;

        //- Return the product mixture
        const thermoType& products(const scalarFieldListSlice&) const;

        //- Reset the mixture to an unburnt state and update EGR
        void reset(PtrList<volScalarField>& Y) const;

        //- Read dictionary
        void read(const dictionary&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inhomogeneousEGRMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "inhomogeneousEGRMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
