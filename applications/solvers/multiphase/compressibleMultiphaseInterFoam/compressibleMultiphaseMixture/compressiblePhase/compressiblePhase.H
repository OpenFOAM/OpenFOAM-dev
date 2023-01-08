/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::compressiblePhase

Description
    Single incompressible phase derived from the phase-fraction.
    Used as part of the multiPhaseMixture for interface-capturing multi-phase
    simulations.

SourceFiles
    compressiblePhase.C

\*---------------------------------------------------------------------------*/

#ifndef compressiblePhase_H
#define compressiblePhase_H

#include "phase.H"
#include "dictionaryEntry.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class compressiblePhase Declaration
\*---------------------------------------------------------------------------*/

class compressiblePhase
:
    public phase
{
    // Private Data

        autoPtr<rhoThermo> thermo_;
        volScalarField Alpha_;
        volScalarField::Internal dgdt_;


public:

    // Constructors

        //- Construct from components
        compressiblePhase
        (
            const word& phaseName,
            const volScalarField& T
        );

        //- Return clone
        autoPtr<compressiblePhase> clone() const;

        //- Return a pointer to a new compressiblePhase created on freestore
        //  from Istream
        class iNew
        {
            const volScalarField& T_;

        public:

            iNew
            (
                const volScalarField& T
            )
            :
                T_(T)
            {}

            autoPtr<compressiblePhase> operator()(Istream& is) const
            {
                return autoPtr<compressiblePhase>
                (
                    new compressiblePhase(is, T_)
                );
            }
        };


    // Member Functions

        //- Return const-access to phase rhoThermo
        const rhoThermo& thermo() const
        {
            return thermo_();
        }

        //- Return access to phase rhoThermo
        rhoThermo& thermo()
        {
            return thermo_();
        }

        //- Return const-access to phase mass-fraction
        const volScalarField& Alpha() const
        {
            return Alpha_;
        }

        //- Return access to phase mass-fraction
        volScalarField& Alpha()
        {
            return Alpha_;
        }

        //- Return const-access to phase divergence
        const volScalarField::Internal& dgdt() const
        {
            return dgdt_;
        }

        //- Return access to phase divergence
        volScalarField::Internal& dgdt()
        {
            return dgdt_;
        }

        void correct
        (
            const volScalarField& p,
            const volScalarField& T
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
