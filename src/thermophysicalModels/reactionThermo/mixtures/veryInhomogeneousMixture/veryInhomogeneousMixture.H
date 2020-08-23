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
    Foam::veryInhomogeneousMixture

Description
    Foam::veryInhomogeneousMixture

SourceFiles
    veryInhomogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef veryInhomogeneousMixture_H
#define veryInhomogeneousMixture_H

#include "basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class veryInhomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class veryInhomogeneousMixture
:
    public basicCombustionMixture
{

public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;

    //- The type of mixture
    typedef ThermoType mixtureType;


private:

    // Private Data

        static const int nSpecies_ = 3;
        static const char* specieNames_[3];

        dimensionedScalar stoicRatio_;

        mixtureType fuel_;
        mixtureType oxidant_;
        mixtureType products_;

        mutable mixtureType mixture_;

        //- Mixture fraction
        volScalarField& ft_;

        //- Fuel mass fraction
        volScalarField& fu_;

        //- Regress variable
        volScalarField& b_;


public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        veryInhomogeneousMixture
        (
            const dictionary&,
            const fvMesh&,
            const word&
        );

        //- Disallow default bitwise copy construction
        veryInhomogeneousMixture
        (
            const veryInhomogeneousMixture<ThermoType>&
        ) = delete;

    //- Destructor
    virtual ~veryInhomogeneousMixture()
    {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "veryInhomogeneousMixture<" + ThermoType::typeName() + '>';
        }

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        const mixtureType& mixture(const scalar, const scalar) const;

        const mixtureType& cellMixture(const label celli) const
        {
            return mixture(ft_[celli], fu_[celli]);
        }

        const mixtureType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                fu_.boundaryField()[patchi][facei]
            );
        }

        const mixtureType& cellReactants(const label celli) const
        {
            return mixture(ft_[celli], ft_[celli]);
        }

        const mixtureType& patchFaceReactants
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture
            (
                ft_.boundaryField()[patchi][facei],
                ft_.boundaryField()[patchi][facei]
            );
        }

        const mixtureType& cellProducts(const label celli) const
        {
            scalar ft = ft_[celli];
            return mixture(ft, fres(ft, stoicRatio().value()));
        }

        const mixtureType& patchFaceProducts
        (
            const label patchi,
            const label facei
        ) const
        {
            scalar ft = ft_.boundaryField()[patchi][facei];
            return mixture(ft, fres(ft, stoicRatio().value()));
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        const ThermoType& specieThermo(const label speciei) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const veryInhomogeneousMixture<ThermoType>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "veryInhomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
