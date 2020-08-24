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
    Foam::inhomogeneousMixture

Description
    Foam::inhomogeneousMixture

SourceFiles
    inhomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef inhomogeneousMixture_H
#define inhomogeneousMixture_H

#include "basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class inhomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class inhomogeneousMixture
:
    public basicCombustionMixture
{

public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;

    //- Mixing type for thermodynamic properties
    typedef ThermoType thermoMixtureType;

    //- Mixing type for transport properties
    typedef ThermoType transportMixtureType;


private:

    // Private Data

        static const int nSpecies_ = 2;
        static const char* specieNames_[2];

        dimensionedScalar stoicRatio_;

        thermoType fuel_;
        thermoType oxidant_;
        thermoType products_;

        mutable thermoType mixture_;

        //- Mixture fraction
        volScalarField& ft_;

        //- Regress variable
        volScalarField& b_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        inhomogeneousMixture(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        inhomogeneousMixture(const inhomogeneousMixture<ThermoType>&) = delete;


    //- Destructor
    virtual ~inhomogeneousMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "inhomogeneousMixture<" + ThermoType::typeName() + '>';
        }

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        const thermoType& mixture(const scalar, const scalar) const;

        const thermoMixtureType& cellThermoMixture(const label celli) const
        {
            return mixture(ft_[celli], b_[celli]);
        }

        const thermoMixtureType& patchFaceThermoMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                b_.boundaryField()[patchi][facei]
            );
        }

        const transportMixtureType& cellTransportMixture
        (
            const label celli
        ) const
        {
            return cellThermoMixture(celli);
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return patchFaceThermoMixture(patchi, facei);
        }

        const transportMixtureType& cellTransportMixture
        (
            const label,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return thermoMixture;
        }

        const transportMixtureType& patchFaceTransportMixture
        (
            const label,
            const label,
            const thermoMixtureType& thermoMixture
        ) const
        {
            return thermoMixture;
        }

        const thermoType& cellReactants(const label celli) const
        {
            return mixture(ft_[celli], 1);
        }

        const thermoType& patchFaceReactants
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                1
            );
        }

        const thermoType& cellProducts(const label celli) const
        {
            return mixture(ft_[celli], 0);
        }

        const thermoType& patchFaceProducts
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                0
            );
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        const ThermoType& specieThermo(const label speciei) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const inhomogeneousMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "inhomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
