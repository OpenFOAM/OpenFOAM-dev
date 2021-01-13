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
    Foam::laminarThermophysicalTransportModels::FickianFourier

Description
    Multi-component Fickian and Fourier based temperature gradient heat flux
    models with optional Soret thermal diffusion of species for laminar flow.

    The mixture diffusion coefficients are specified as Function2<scalar>s of
    pressure and temperature but independent of composition.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

Usage
    \verbatim
    laminar
    {
        model           FickianFourier;

        mixtureDiffusionCoefficients yes;

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
    laminar
    {
        model           FickianFourier;

        mixtureDiffusionCoefficients no;

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
    FickianFourier.C

\*---------------------------------------------------------------------------*/

#include "Fickian.H"
#include "unityLewisFourier.H"

#ifndef FickianFourier_H
#define FickianFourier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class FickianFourier Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class FickianFourier
:
    public Fickian
    <
        unityLewisFourier<laminarThermophysicalTransportModel>
    >
{

protected:

    // Protected data

        // Model coefficients

            //- Turbulent Schmidt number []
            dimensionedScalar Sct_;


public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("FickianFourier");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        FickianFourier
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~FickianFourier()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FickianFourier.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
