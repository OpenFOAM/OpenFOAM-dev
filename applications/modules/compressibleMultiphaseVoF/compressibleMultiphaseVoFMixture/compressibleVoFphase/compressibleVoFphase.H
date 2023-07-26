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
    Foam::compressibleVoFphase

Description
    Single compressible phase derived from the VoFphase.

    Used in compressibleMultiphaseVoFMixture for multiphase
    interface-capturing simulations.

SourceFiles
    compressibleVoFphase.C

See also
    Foam::VoFphase
    Foam::compressibleMultiphaseVoFMixture

\*---------------------------------------------------------------------------*/

#ifndef compressibleVoFphase_H
#define compressibleVoFphase_H

#include "VoFphase.H"
#include "rhoFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class compressibleVoFphase Declaration
\*---------------------------------------------------------------------------*/

class compressibleVoFphase
:
    public VoFphase
{
    // Private Data

        //- Phase thermo
        autoPtr<rhoFluidThermo> thermo_;

        //- Phase mass-fraction
        volScalarField Alpha_;

        //- Phase compressibility contribution
        volScalarField::Internal dgdt_;


public:

    // Constructors

        //- Construct from components
        compressibleVoFphase
        (
            const word& name,
            const fvMesh& mesh,
            const volScalarField& T
        );

        //- Return clone
        virtual autoPtr<VoFphase> clone() const;

        //- Return a pointer to a new compressibleVoFphase
        //  created on freestore from Istream
        class iNew
        :
            public VoFphase::iNew
        {
            const volScalarField& T_;

        public:

            iNew(const fvMesh& mesh, const volScalarField& T)
            :
                VoFphase::iNew(mesh),
                T_(T)
            {}

            virtual autoPtr<VoFphase> operator()(Istream& is) const
            {
                const word name(is);

                return autoPtr<VoFphase>
                (
                    new compressibleVoFphase(name, mesh_, T_)
                );
            }
        };


    // Member Functions

        //- Return const-access to phase rhoFluidThermo
        const rhoFluidThermo& thermo() const
        {
            return thermo_();
        }

        //- Return access to phase rhoFluidThermo
        rhoFluidThermo& thermo()
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
