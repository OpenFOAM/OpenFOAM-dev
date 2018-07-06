/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::combustionModels::FSD

Description

    Flame Surface Dennsity (FDS) combustion model.

    The fuel source term is given by mgft*pc*omegaFuelBar.

    where:
          mgft: filtered flame area.
          pc:   probability of the combustion progress.
          omegaFuelBar: filtered consumption speed per unit of flame area.

    pc is considered from the IFC solution.
    omegaFuelBar is calculated solving a relaxation equation which tends to
    omegaEq. This omegaEq is obtained from the flamelet solution for
    different strain rates and fit using a expential distribution.

    The spacial distribution of the consumption speed (omega) is obtained also
    from a strained flamelet solution and it is assumed to have a gaussian
    distribution.

    If the grid resolution is not enough to resolve the flame, the consumption
    speed distribution is linearly thickened conserving the overall heat
    release.

    If the turbulent fluctuation of the mixture fraction at the sub-grid level
    is large (>1e-04) then a beta pdf is used for filtering.

    At the moment the flame area combustion model is only fit to work in a LES
    frame work. In RAS the subgrid fluctuation has to be solved by an extra
    transport equation.

SourceFiles
    FSD.C

\*---------------------------------------------------------------------------*/

#ifndef FSD_H
#define FSD_H

#include "singleStepCombustion.H"
#include "reactionRateFlameArea.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                          Class FSD Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class FSD
:
    public singleStepCombustion <ReactionThermo, ThermoType>
{
    // Private data

        //- Auto pointer to consumption speed per unit of flame area model
        autoPtr<reactionRateFlameArea> reactionRateFlameArea_;

        //- Mixture fraction
        volScalarField ft_;

        //- Fuel mass concentration on the fuel stream
        dimensionedScalar YFuelFuelStream_;

        //- Oxygen mass concentration on the oxydizer stream
        dimensionedScalar YO2OxiStream_;

        //- Similarity constant for the sub-grid ft fluctuations
        scalar Cv_;

        //- Model constant
        scalar C_;

        //- Lower flammability limit
        scalar ftMin_;

        //- Upper flammability limit
        scalar ftMax_;

        //- Dimension of the ft space. Used to integrate the beta-pdf
        scalar ftDim_;

        //- Minimum mixture freaction variance to calculate pdf
        scalar ftVarMin_;


    // Private Member Functions

        //- Calculate the normalised fuel source term
        void calculateSourceNorm();

        //- Disallow copy construct
        FSD(const FSD&);

        //- Disallow default bitwise assignment
        void operator=(const FSD&);


public:

    //- Runtime type information
    TypeName("FSD");


    // Constructors

        //- Construct from components
        FSD
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~FSD();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Update properties
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FSD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
