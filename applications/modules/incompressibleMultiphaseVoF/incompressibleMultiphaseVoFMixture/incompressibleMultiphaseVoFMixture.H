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
    Foam::incompressibleMultiphaseVoFMixture

Description
    Incompressible multiphase mixture for interface-capturing simulations.

    Derived from viscosity so that it can be used in conjunction with the
    incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface between each
    phase-pair.

SourceFiles
    incompressibleMultiphaseVoFMixture.C

See also
    Foam::multiphaseVoFMixture
    Foam::incompressibleVoFphase

\*---------------------------------------------------------------------------*/

#ifndef incompressibleMultiphaseVoFMixture_H
#define incompressibleMultiphaseVoFMixture_H

#include "multiphaseVoFMixture.H"
#include "viscosity.H"
#include "UPtrListDictionary.H"
#include "incompressibleVoFphase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompressibleMultiphaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleMultiphaseVoFMixture
:
    public multiphaseVoFMixture,
    public viscosity
{
    // Private Data

        //- Dictionary of phases
        UPtrListDictionary<incompressibleVoFphase> phases_;

        //- Mixture density
        volScalarField rho_;

        //- Mixture viscosity
        volScalarField nu_;


    // Private Member Functions

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;


public:

    TypeName("incompressibleMultiphaseVoFMixture");


    // Constructors

        //- Construct from components
        incompressibleMultiphaseVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~incompressibleMultiphaseVoFMixture()
    {}


    // Member Functions

        //- Return the phases
        UPtrListDictionary<incompressibleVoFphase>& phases()
        {
            return phases_;
        }

        //- Return the mixture density
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Correct the mixture properties
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
