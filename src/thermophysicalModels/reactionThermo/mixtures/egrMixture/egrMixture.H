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
    Foam::egrMixture

Description
    Foam::egrMixture

SourceFiles
    egrMixture.C

\*---------------------------------------------------------------------------*/

#ifndef egrMixture_H
#define egrMixture_H

#include "basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class egrMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class egrMixture
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

        static const int nSpecies_ = 3;
        static const char* specieNames_[3];

        dimensionedScalar stoicRatio_;

        ThermoType fuel_;
        ThermoType oxidant_;
        ThermoType products_;

        mutable thermoType mixture_;

        //- Mixture fraction
        volScalarField& ft_;

        //- Regress variable
        volScalarField& b_;

        //- Residual gases
        volScalarField& egr_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phaseName
        egrMixture(const dictionary&, const fvMesh&, const word&);

        //- Disallow default bitwise copy construction
        egrMixture(const egrMixture<ThermoType>&) = delete;


    //- Destructor
    virtual ~egrMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "egrMixture<" + ThermoType::typeName() + '>';
        }

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        const thermoType& mixture
        (
            const scalar,
            const scalar,
            const scalar
        ) const;

        const thermoMixtureType& cellThermoMixture(const label celli) const
        {
            return mixture(ft_[celli], b_[celli], egr_[celli]);
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
                b_.boundaryField()[patchi][facei],
                egr_.boundaryField()[patchi][facei]
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
            return mixture(ft_[celli], 1, egr_[celli]);
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
                1,
                egr_.boundaryField()[patchi][facei]
            );
        }

        const thermoType& cellProducts(const label celli) const
        {
            return mixture(ft_[celli], 0, 0);
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
                0,
                0
            );
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        const ThermoType& specieThermo(const label speciei) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const egrMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "egrMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
