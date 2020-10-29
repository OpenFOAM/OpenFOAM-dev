/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::turbulenceThermophysicalTransportModels::nonUnityLewisEddyDiffusivity

Description
    Non-unity-Lewis-Eddy-diffusivity based temperature gradient heat flux model
    for RAS or LES of turbulent flow. Allows independent specification of
    turbulent Prandtl and Schmidt numbers.  Unity laminar Lewis number is
    assumed.

    The heat flux source is implemented as an implicit energy correction to the
    temperature gradient based flux source.  At convergence the energy
    correction is 0.

Usage
    \verbatim
    LES
    {
        model           nonUnityLewisEddyDiffusivity;
        Prt             0.85;
        Sct             0.7;
    }
    \endverbatim

SourceFiles
    nonUnityLewisEddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#include "unityLewisEddyDiffusivity.H"

#ifndef nonUnityLewisEddyDiffusivity_H
#define nonUnityLewisEddyDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class nonUnityLewisEddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceThermophysicalTransportModel>
class nonUnityLewisEddyDiffusivity
:
    public unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
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
    TypeName("nonUnityLewisEddyDiffusivity");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        nonUnityLewisEddyDiffusivity
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~nonUnityLewisEddyDiffusivity()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusivity for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const
        {
            return volScalarField::New
            (
                "DEff",
                this->thermo().alphaEff((this->Prt_/Sct_)*this->alphat())
            );
        }

        //- Effective mass diffusivity for a given specie mass-fraction
        //  for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const
        {
            return this->thermo().alphaEff
            (
                this->Prt_.value()/Sct_.value()*this->alphat(patchi),
                patchi
            );
        }

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonUnityLewisEddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
