/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Nucleation source due to reactions. Applicable for irreversible
    reactions only.

SourceFiles
    reactionDriven.C

\*---------------------------------------------------------------------------*/

#ifndef reactionDriven_H
#define reactionDriven_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
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

        //- Sphere-equivalent nucleation diameter
        const dimensionedScalar dNuc_;

        //- Phase in which the reaction occurs
        const phaseModel& reactingPhase_;

        //- Name of the mass transfer rate
        const word dmdtfName_;

        //- Name of the specie which nucleates
        const word specieName_;


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

        //- Add to nucleationRate
        virtual void addToNucleationRate
        (
            volScalarField& nucleationRate,
            const label i
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
