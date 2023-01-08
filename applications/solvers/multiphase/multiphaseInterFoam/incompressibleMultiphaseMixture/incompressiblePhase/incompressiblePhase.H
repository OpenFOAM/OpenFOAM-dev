/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::incompressiblePhase

Description
    Single incompressible phase derived from the phase.

    Used as part of the multiphaseMixture for interface-capturing
    multiphase simulations.

SourceFiles
    incompressiblePhase.C

\*---------------------------------------------------------------------------*/

#ifndef incompressiblePhase_H
#define incompressiblePhase_H

#include "phase.H"
#include "dictionaryEntry.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class incompressiblePhase Declaration
\*---------------------------------------------------------------------------*/

class incompressiblePhase
:
    public phase
{
    // Private Data

        autoPtr<viscosityModel> nuModel_;
        dimensionedScalar rho_;


public:

    // Constructors

        //- Construct from components
        incompressiblePhase
        (
            const word& name,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<incompressiblePhase> clone() const;

        //- Return a pointer to a new incompressiblePhase created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<incompressiblePhase> operator()(Istream& is) const
            {
                const word name(is);

                return autoPtr<incompressiblePhase>
                (
                    new incompressiblePhase(name, mesh_)
                );
            }
        };


    // Member Functions

        //- Return const-access to incompressiblePhase1 viscosityModel
        const viscosityModel& nuModel() const
        {
            return nuModel_();
        }

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nuModel_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nuModel_->nu(patchi);
        }

        //- Return const-access to incompressiblePhase1 density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }

        //- Correct the incompressiblePhase properties
        void correct();

        //- Inherit read from phase
        using phase::read;

        //- Read base incompressiblePhaseProperties dictionary
        bool read(const dictionary& incompressiblePhaseDict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
