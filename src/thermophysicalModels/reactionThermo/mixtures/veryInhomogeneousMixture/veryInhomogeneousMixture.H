/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class veryInhomogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class veryInhomogeneousMixture
:
    public basicMultiComponentMixture
{
    // Private data

        static const int nSpecies_ = 3;
        static const char* specieNames_[3];

        dimensionedScalar stoicRatio_;

        ThermoType fuel_;
        ThermoType oxidant_;
        ThermoType products_;

        mutable ThermoType mixture_;

        //- Mixture fraction
        volScalarField& ft_;

        //- Fuel mass fraction
        volScalarField& fu_;

        //- Regress variable
        volScalarField& b_;

        //- Construct as copy (not implemented)
        veryInhomogeneousMixture(const veryInhomogeneousMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary and mesh
        veryInhomogeneousMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~veryInhomogeneousMixture()
    {}


    // Member functions

        const dimensionedScalar& stoicRatio() const
        {
            return stoicRatio_;
        }

        const ThermoType& mixture(const scalar, const scalar) const;

        const ThermoType& cellMixture(const label celli) const
        {
            return mixture(ft_[celli], fu_[celli]);
        }

        const ThermoType& patchFaceMixture
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

        const ThermoType& cellReactants(const label celli) const
        {
            return mixture(ft_[celli], ft_[celli]);
        }

        const ThermoType& patchFaceReactants
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

        const ThermoType& cellProducts(const label celli) const
        {
            scalar ft = ft_[celli];
            return mixture(ft, fres(ft, stoicRatio().value()));
        }

        const ThermoType& patchFaceProducts
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
        const ThermoType& getLocalThermo(const label specieI) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "veryInhomogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
