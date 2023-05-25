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
    Foam::compressibleMultiphaseVoFMixture

Description
    Compressible multiphase mixture for interface-capturing simulations.

    Derived from viscosity so that it can be used in conjunction with the
    compressible turbulence models.

    Surface tension and contact-angle is handled for the interface between each
    phase-pair.

SourceFiles
    compressibleMultiphaseVoFMixture.C

See also
    Foam::multiphaseVoFMixture
    Foam::compressibleVoFphase

\*---------------------------------------------------------------------------*/

#ifndef compressibleMultiphaseVoFMixture_H
#define compressibleMultiphaseVoFMixture_H

#include "compressibleMultiphaseVoFMixtureThermo.H"
#include "multiphaseVoFMixture.H"
#include "viscosity.H"
#include "UPtrListDictionary.H"
#include "compressibleVoFphase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class compressibleMultiphaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class compressibleMultiphaseVoFMixture
:
    public compressibleMultiphaseVoFMixtureThermo,
    public multiphaseVoFMixture,
    public viscosity
{
    // Private Data

        //- Dictionary of phases
        UPtrListDictionary<compressibleVoFphase> phases_;

        volScalarField rho_;


public:

    TypeName("compressibleMultiphaseVoFMixture");


    // Constructors

        //- Construct from fvMesh
        compressibleMultiphaseVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~compressibleMultiphaseVoFMixture()
    {}


    // Member Functions

        //- Return the phases
        UPtrListDictionary<compressibleVoFphase>& phases()
        {
            return phases_;
        }

        //- Return the mixture density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return true if all phases are incompressible
        bool incompressible() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Return the mixture compressibility/density
        virtual tmp<volScalarField> psiByRho() const;

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();

        //- Return the effective temperature transport coefficient
        //  derived from the phase internal energy equations i.e. <kappa/Cv>
        virtual tmp<volScalarField> alphaEff
        (
            const volScalarField& nut
        ) const;

        //- Return the phase-averaged reciprocal Cv
        tmp<volScalarField> rCv() const;

        //- Update properties
        virtual void correct();

        //- Update densities for given pressure change
        void correctRho(const volScalarField& dp);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
