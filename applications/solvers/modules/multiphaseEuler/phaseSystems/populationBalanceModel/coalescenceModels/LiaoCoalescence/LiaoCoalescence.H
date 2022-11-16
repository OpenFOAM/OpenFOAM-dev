/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::diameterModels::coalescenceModels::LiaoCoalescence

Description
    Bubble coalescence model of Liao et al. (2015). The terminal velocities and
    drag coefficients are computed by an iterative procedure based on the drag
    model of Ishii and Zuber (1979) at the beginning of the simulation, assuming
    single bubbles rising in quiescent liquid.

    Note that the original article contains a mistake concerning the value of
    the coefficient CEff. A value of 2.5 instead of 5.0 should be used when
    following the Weber number definition in the paper.

    References:
    \verbatim
        Liao, Y.; Rzehak, R.; Lucas, D.; Krepper, E. (2015).
        Baseline closure models for dispersed bubbly flow: Bubble coalescence
        and breakup.
        Chemical Engineering Science, 122, 336-349.

        Liao, Y., Rzehak, R., Lucas, D., & Krepper, E. (2021).
        Corrigendum to "Baseline closure model for dispersed bubbly flow:
        Bubble coalescence and breakup" [Chem. Eng. Sci. 122 (2015) 336–349].
        Chemical Engineering Science, 241, 116708.

        Ishii, M., & Zuber, N. (1979).
        Drag coefficient and relative velocity in bubbly, droplet or particulate
        flows.
        AIChE Journal, 25(5), 843-855.
    \endverbatim

Usage
    \table
        Property     | Description                 | Required  | Default value
        PMax         | Maximum packing limit       | no        | 0.8
        CPackMax     | Maximum CPack coefficient   | no        | 1.0e5
        AH           | Hamaker constant            | no        | 3.7e-20
        CEff         | coefficient CEff            | no        | 2.5
        CTurb        | coefficient CTurb           | no        | 1.0
        CBuoy        | coefficient CBuoy           | no        | 1.0
        CShear       | coefficient CShear          | no        | 1.0
        CEddy        | coefficient CEddy           | no        | 1.0
        CWake        | coefficient CWake           | no        | 1.0
        turbulence   | Switch for turbulence       | yes       | none
        buoyancy     | Switch for buoyancy         | yes       | none
        laminarShear | Switch for laminar shear    | yes       | none
        eddyCapture  | Switch for eddy capture     | yes       | none
        wake         | Switch for wake entrainment | yes       | none
    \endtable

SourceFiles
    LiaoCoalescence.C

\*---------------------------------------------------------------------------*/

#ifndef LiaoCoalescence_H
#define LiaoCoalescence_H

#include "coalescenceModel.H"
#include "LiaoBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace coalescenceModels
{

/*---------------------------------------------------------------------------*\
                       Class LiaoCoalescence Declaration
\*---------------------------------------------------------------------------*/

class LiaoCoalescence
:
    public coalescenceModel,
    public LiaoBase
{
    // Private Data

        //- Maximum packing limit of bubbles
        dimensionedScalar PMax_;

        //- Hamaker constant AH
        dimensionedScalar AH_;

        //- Optional coefficient CEff
        dimensionedScalar CEff_;

        //- Optional coefficient CTurb
        dimensionedScalar CTurb_;

        //- Optional coefficient CBuoy
        dimensionedScalar CBuoy_;

        //- Optional coefficient CShear
        dimensionedScalar CShear_;

        //- Optional coefficient CEddy
        dimensionedScalar CEddy_;

        //- Optional coefficient CWake
        dimensionedScalar CWake_;

        //- Switch for turbulent collisions
        Switch turbulence_;

        //- Switch for buoyancy-induced collisions
        Switch buoyancy_;

        //- Switch for velocity gradient-induced collisions
        Switch laminarShear_;

        //- Switch for eddy capture-induced collisions
        Switch eddyCapture_;

        //- Switch for wake entrainment-induced collisions
        Switch wakeEntrainment_;

        //- Coefficient containing maximum packing limit
        volScalarField CPack_;

        //- Limiter for the CPack coefficient
        dimensionedScalar CPackMax_;

        //- Critical diameter
        volScalarField dCrit_;

        //- Relative velocity for turbulent collisions
        volScalarField uRelTurb_;

        //- Relative velocity for buoyancy-induced collisions
        volScalarField uRelBuoy_;

        //- Relative velocity for shear-induced collisions
        volScalarField uRelShear_;


public:

    //- Runtime type information
    TypeName("Liao");

    // Constructor

        LiaoCoalescence
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~LiaoCoalescence()
    {}


    // Member Functions

        //- Precompute diameter independent expressions
        virtual void precompute();

        //- Add to coalescenceRate
        virtual void addToCoalescenceRate
        (
            volScalarField& coalescenceRate,
            const label i,
            const label j
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
