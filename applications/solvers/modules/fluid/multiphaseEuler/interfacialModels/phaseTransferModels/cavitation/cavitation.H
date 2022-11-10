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
    Foam::phaseTransferModels::cavitation

Description
    A phase transfer model which represents cavitation.

    Example usage:
    \verbatim
    phaseTransfer
    {
        vapour_dispersedIn_water
        {
            type        cavitation;

            model       Kunz;

            liquid      water;

            pSat        80000;

            UInf        5.33;
            tInf        0.028142589;
            Cc          100;
            Cv          100;
        }
    }
    \endverbatim

SourceFiles
    cavitation.C

\*---------------------------------------------------------------------------*/

#ifndef cavitation_H
#define cavitation_H

#include "phaseTransferModel.H"
#include "cavitationModel.H"
#include "phaseModel.H"
#include "dispersedPhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace phaseTransferModels
{

/*---------------------------------------------------------------------------*\
                       Class cavitation Declaration
\*---------------------------------------------------------------------------*/

class cavitation
:
    public phaseTransferModel
{
private:

    // Private Data

        //- Interface
        const phaseInterface interface_;

        //- Cavitation model
        autoPtr<compressible::cavitationModel> cavitation_;


public:

    //- Runtime type information
    TypeName("cavitation");


    // Constructors

        //- Construct from components
        cavitation
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~cavitation();


    // Member Functions

        //- Whether or not there is a mass transfer rate for the mixture
        virtual bool mixture() const;

        //- The mass transfer rate for the mixture
        virtual tmp<volScalarField> dmdtf() const;

        //- The derivative of mass transfer w.r.t. pressure
        virtual tmp<volScalarField> d2mdtdpf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
