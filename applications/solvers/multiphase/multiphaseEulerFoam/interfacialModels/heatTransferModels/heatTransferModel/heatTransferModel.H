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
    Foam::heatTransferModel

Description

SourceFiles
    heatTransferModel.C
    newHeatTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferModel_H
#define heatTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class heatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferModel
{
protected:

    // Protected data

        const dictionary& interfaceDict_;
        const volScalarField& alpha1_;
        const phaseModel& phase1_;
        const phaseModel& phase2_;


public:

    //- Runtime type information
    TypeName("heatTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferModel,
            dictionary,
            (
                const dictionary& interfaceDict,
                const volScalarField& alpha1,
                const phaseModel& phase1,
                const phaseModel& phase2
            ),
            (interfaceDict, alpha1, phase1, phase2)
        );


    // Constructors

        heatTransferModel
        (
            const dictionary& interfaceDict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    virtual ~heatTransferModel();


    // Selectors

        static autoPtr<heatTransferModel> New
        (
            const dictionary& interfaceDict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    // Member Functions

        //- The heat-transfer function K used in the enthalpy eq.
        //    ddt(alpha1*rho1*ha) + ... = ... alpha1*alpha2*K*(Ta - Tb)
        //    ddt(alpha2*rho2*hb) + ...  = ... alpha1*alpha2*K*(Tb - Ta)
        // ********************************** NB! *****************************
        // for numerical reasons alpha1 and alpha2 has been
        // extracted from the heat-transfer function K,
        // so you MUST divide K by alpha1*alpha2 when implementing the
        // heat-transfer function
        // ********************************** NB! *****************************
        virtual tmp<volScalarField> K(const volScalarField& Ur) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
