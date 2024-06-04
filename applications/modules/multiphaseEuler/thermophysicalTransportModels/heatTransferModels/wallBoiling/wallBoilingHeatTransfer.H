/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
    Foam::heatTransferModels::wallBoiling

Description
    A heat transfer model for simulation of sub-cooled nucleate wall boiling on
    the surface of a third phase with runtime selectable sub-models. This model
    should be applied to the liquid side of the liquid-third-phase interface.

Usage
    Settings below are very similar to that of the
    alphatWallBoilingWallFunction (refer to this model's documentation for
    further details). The only differences are that an underlying heat transfer
    model is also required, and the phase-identifying entries are somewhat
    different; the liquid and vapour phases must be specified explicitly.

    A corresponding heat transfer model on the other side of the
    liquid-third-phase interface is required by this model, and is also shown
    below.

    Example usage:
    \verbatim
    bed_dispersedIn_liquid_inThe_liquid
    {
        type            wallBoiling;

        vapourPhase     gas;

        heatTransferModel
        {
            type            Gunn;
        }

        partitioningModel
        {
            type            Lavieville; // phaseFraction, linear, cosine
            alphaCrit       0.2;
        }
        nucleationSiteModel
        {
            type            LemmertChawla; // KocamustafaogullariIshii
        }
        departureDiameterModel
        {
            type            TolubinskiKostanchuk; // KocamustafaogullariIshii
        }
        departureFrequencyModel
        {
            type            KocamustafaogullariIshii; // Cole
            Cf              1.18;
        }
    }

    bed_dispersedIn_liquid_inThe_bed
    {
        type            spherical;
    }
    \endverbatim

See also
    Foam::compressible::alphatWallBoilingWallFunctionFvPatchScalarField

SourceFiles
    wallBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoilingHeatTransfer_H
#define wallBoilingHeatTransfer_H

#include "heatTransferModel.H"
#include "partitioningModel.H"
#include "nucleationSiteModel.H"
#include "departureDiameterModel.H"
#include "departureFrequencyModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                   Class wallBoilingHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class wallBoilingHeatTransfer
:
    public heatTransferModel
{
    // Private Data

        //- Interface
        const dispersedPhaseInterface interface_;

        //- Interface on the other side
        autoPtr<phaseInterface> otherInterface_;

        //- Name of the vapour phase
        const word vapourPhaseName_;

        //- Pointer to the underlying heat transfer model
        autoPtr<heatTransferModel> heatTransferModel_;

        //- Relaxation factor
        const scalar relax_;

        //- Run-time selected heat flux partitioning model
        autoPtr<Foam::wallBoilingModels::partitioningModel>
            partitioningModel_;

        //- Run-time selected nucleation site density model
        autoPtr<Foam::wallBoilingModels::nucleationSiteModel>
            nucleationSiteModel_;

        //- Run-time selected bubble departure diameter model
        autoPtr<Foam::wallBoilingModels::departureDiameterModel>
            departureDiamModel_;

        //- Run-time selected bubble departure frequency model
        autoPtr<Foam::wallBoilingModels::departureFrequencyModel>
            departureFreqModel_;

        //- At wall liquid phase fraction
        mutable volScalarField wetFraction_;

        //- Bubble departure diameter
        mutable volScalarField dDep_;

        //- Bubble departre frequency
        mutable volScalarField fDep_;

        //- Nucleation site density
        mutable volScalarField nucleationSiteDensity_;

        //- Evapouration mass transfer rate
        mutable volScalarField dmdtf_;

        //- Quenching heat transfer rate
        mutable volScalarField qq_;

        //- Wall Temperature
        mutable volScalarField Tsurface_;

        //- Heat transfer coefficient (needed because lagging)
        mutable volScalarField K_;


public:

    //- Runtime type information
    TypeName("wallBoiling");


    // Constructors

        //- Construct from components
        wallBoilingHeatTransfer
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~wallBoilingHeatTransfer();


    // Member Functions

        //- The heat transfer function K used in the enthalpy equation
        tmp<volScalarField> K(const scalar residualAlpha) const;

        //- Is there phase change mass transfer for this phaseInterface
        bool activePhaseInterface(const phaseInterfaceKey&) const;

        //- True if the sign of dmdtf should be changed
        bool flipSign() const;

        //- Return the rate of phase-change
        const volScalarField& dmdtf() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
