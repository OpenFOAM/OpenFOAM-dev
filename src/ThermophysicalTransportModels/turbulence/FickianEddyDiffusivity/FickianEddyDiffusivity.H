/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::turbulenceThermophysicalTransportModels::FickianEddyDiffusivity

Description
    Multi-component Fickian and eddy-diffusivity turbulent based temperature
    gradient heat flux model for RAS or LES of turbulent flow with optional
    Soret thermal diffusion of species.

    The mixture diffusion coefficients are specified as Function2<scalar>s of
    pressure and temperature but independent of composition.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

Usage
    \verbatim
    RAS
    {
        model           FickianEddyDiffusivity;

        mixtureDiffusionCoefficients yes;

        Prt             0.85;
        Sct             0.7;

        Dm // [m^2/s]
        {
            O2 1e-2;
            O3 5e-2;
            N2 1e-2;
        }

        DT // [kg/m/s] Optional
        {
            O2 1e-2;
            O3 5e-2;
            N2 1e-2;
        }
    }
    \endverbatim

    or if binary mass diffusion coefficient functions are available they can be
    mixed to form the mass diffusion coefficients w.r.t. the mixture:

    \verbatim
    RAS
    {
        model           FickianEddyDiffusivity;

        mixtureDiffusionCoefficients no;

        Prt             0.85;
        Sct             0.7;

        D // [m^2/s]
        {
            O2-O2 1e-2;
            O3-O3 5e-2;
            N2-N2 1e-2;
            O3-O2 5e-2;
            O3-N2 5e-2;
            O2-N2 1e-2;
        }

        DT // [kg/m/s] Optional
        {
            O2 1e-2;
            O3 5e-2;
            N2 1e-2;
        }
    }
    \endverbatim

SourceFiles
    FickianEddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#include "Fickian.H"
#include "unityLewisEddyDiffusivity.H"

#ifndef FickianEddyDiffusivity_H
#define FickianEddyDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class FickianEddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceThermophysicalTransportModel>
class FickianEddyDiffusivity
:
    public Fickian
    <
        unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
    >
{

protected:

    // Protected data

        // Model coefficients

            //- Turbulent Schmidt number []
            dimensionedScalar Sct_;


public:

    typedef typename TurbulenceThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        TurbulenceThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename TurbulenceThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("FickianEddyDiffusivity");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        FickianEddyDiffusivity
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~FickianEddyDiffusivity()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FickianEddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
