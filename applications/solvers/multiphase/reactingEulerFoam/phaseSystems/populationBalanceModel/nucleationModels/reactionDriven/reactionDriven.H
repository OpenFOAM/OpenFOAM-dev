/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::diameterModels::nucleationModels::constantNucleation

Description
    Applicable for irreversible reactions only.

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef reactionDriven_H
#define reactionDriven_H

#include "nucleationModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "ReactingPhaseModel.H"
#include "MovingPhaseModel.H"
#include "ThermoPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        ReactingPhaseModel
        <
            MovingPhaseModel
            <
                ThermoPhaseModel<phaseModel, rhoReactionThermo>
            >,
            CombustionModel<rhoReactionThermo>
        >
        reactingPhaseModel;

namespace diameterModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                       Class reactionDriven Declaration
\*---------------------------------------------------------------------------*/

class reactionDriven
:
    public nucleationModel
{
    // Private Data

        //- Sphere-eqivalent nucleation diameter
        dimensionedScalar dNuc_;

        //- Velocity group in which the nucleation occurs
        const velocityGroup& velGroup_;

        //- Phase in which the reaction occurs
        const reactingPhaseModel& reactingPhase_;

        //- Product species of the reaction
        volScalarField& specie_;

        //- Nucleation mass transfer rate
        volScalarField nDmdt_;


public:

    //- Runtime type information
    TypeName("reactionDriven");

    // Constructor

        reactionDriven
        (
            const populationBalanceModel& popBal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reactionDriven()
    {}


    // Member Functions

        //- Switch which determines, whether the nucleation model introduce mass
        //  into the domain through the populationBalanceModel class, i.e. pDmdt
        virtual bool pDmDt() const
        {
            return true;
        };

        //- If the nucleation model uses pDmdt, it has to ensure that the
        //  relevant phasePair exists in the speciesDmdt-table
        virtual void registerPair
        (
            populationBalanceModel::speciesDmdtTable& speciesDmdt
        ) const;

        //- Correct diameter independent expressions
        virtual void correct();

        //- Add to nucleationRate
        virtual void addToNucleationRate
        (
            volScalarField& nucleationRate,
            const label i
        );

        //- Add the corresponding species sources
        virtual void addToSpeciesDmDt
        (
            populationBalanceModel::speciesDmdtTable& speciesDmdt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
