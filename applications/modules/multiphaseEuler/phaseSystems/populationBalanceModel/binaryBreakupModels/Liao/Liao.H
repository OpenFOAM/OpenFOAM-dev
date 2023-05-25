/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::diameterModels::binaryBreakupModels::Liao

Description
    Bubble breakup model of Liao et al. (2015). The terminal velocities and drag
    coefficients are computed by an iterative procedure based on the drag model
    of Ishii and Zuber (1979) at the beginning of the simulation, assuming
    single bubbles rising in quiescent liquid.

    References:
    \verbatim
        Liao, Y.; Rzehak, R.; Lucas, D.; Krepper, E. (2015).
        Baseline closure models for dispersed bubbly flow: Bubble coalescence
        and breakup.
        Chemical Engineering Science, 122, 336-349.

        Ishii, M., & Zuber, N. (1979).
        Drag coefficient and relative velocity in bubbly, droplet or particulate
        flows.
        AIChE Journal, 25(5), 843-855.
    \endverbatim

Usage
    \table
        Property       | Description              | Required    | Default value
        BTurb          | coefficient BTurb        | no          | 1.0
        BShear         | coefficient Bbuoy        | no          | 1.0
        BEddy          | coefficient BEddy        | no          | 1.0
        BFric          | coefficient Bwake        | no          | 1.0
        turbulence     | Switch for turbulence    | yes         | none
        laminarShear   | Switch for laminar shear | yes         | none
        turbulentShear | Switch for eddies        | yes         | none
        friction       | Switch for friction      | yes         | none
    \endtable

SourceFiles
    Liao.C

\*---------------------------------------------------------------------------*/

#ifndef Liao_H
#define Liao_H

#include "binaryBreakupModel.H"
#include "LiaoBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace binaryBreakupModels
{

/*---------------------------------------------------------------------------*\
                            Class Liao Declaration
\*---------------------------------------------------------------------------*/

class Liao
:
    public binaryBreakupModel,
    public LiaoBase
{
    // Private Data

        //- Optional coefficient BTurb
        dimensionedScalar BTurb_;

        //- Optional coefficient BShear
        dimensionedScalar BShear_;

        //- Optional coefficient BEddy
        dimensionedScalar BEddy_;

        //- Optional coefficient BFric
        dimensionedScalar BFric_;

        //- Switch for considering turbulent velocity fluctuation
        Switch turbulence_;

        //- Switch for considering velocity gradients in the bulk
        Switch laminarShear_;

        //- Switch for considering eddies
        Switch turbulentShear_;

        //- Switch for considering interfacial friction
        Switch interfacialFriction_;


public:

    //- Runtime type information
    TypeName("Liao");

    // Constructor

        Liao
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Liao()
    {}


    // Member Functions

        //- Precompute diameter independent expressions
        virtual void precompute();

        //- Add to binary breakupRate
        virtual void addToBinaryBreakupRate
        (
            volScalarField& binaryBreakupRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace binaryBreakupModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
