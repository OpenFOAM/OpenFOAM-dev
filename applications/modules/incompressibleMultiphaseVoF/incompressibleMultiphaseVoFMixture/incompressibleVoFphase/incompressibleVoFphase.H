/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::incompressibleVoFphase

Description
    Single incompressible phase derived from the VoFphase.

    Used in incompressibleMultiphaseVoFMixture for multiphase
    interface-capturing simulations.

SourceFiles
    incompressibleVoFphase.C

See also
    Foam::VoFphase
    Foam::incompressibleMultiphaseVoFMixture

\*---------------------------------------------------------------------------*/

#ifndef incompressibleVoFphase_H
#define incompressibleVoFphase_H

#include "VoFphase.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class incompressibleVoFphase Declaration
\*---------------------------------------------------------------------------*/

class incompressibleVoFphase
:
    public VoFphase
{
    // Private Data

        //- Phase viscosity model
        autoPtr<viscosityModel> nuModel_;

        //- Phase density
        dimensionedScalar rho_;


public:

    // Constructors

        //- Construct from components
        incompressibleVoFphase
        (
            const word& name,
            const fvMesh& mesh
        );

        //- Return clone
        virtual autoPtr<VoFphase> clone() const;

        //- Return a pointer to a new incompressibleVoFphase
        //  created on freestore from Istream
        class iNew
        :
            public VoFphase::iNew
        {
        public:

            iNew(const fvMesh& mesh)
            :
                VoFphase::iNew(mesh)
            {}

            virtual autoPtr<VoFphase> operator()(Istream& is) const
            {
                const word name(is);

                return autoPtr<VoFphase>
                (
                    new incompressibleVoFphase(name, mesh_)
                );
            }
        };


    // Member Functions

        //- Return const-access to incompressibleVoFphase1 viscosityModel
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

        //- Return const-access to incompressibleVoFphase1 density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }

        //- Correct the incompressibleVoFphase properties
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
