/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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
    Foam::nonlinearEddyViscosity

Description
    Eddy viscosity turbulence model with non-linear correction base class

SourceFiles
    nonlinearEddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef nonlinearEddyViscosity_H
#define nonlinearEddyViscosity_H

#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonlinearEddyViscosity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class nonlinearEddyViscosity
:
    public eddyViscosity<BasicMomentumTransportModel>
{
    // Private member functions

        //- Return the source term for the momentum equation
        template<class RhoFieldType>
        tmp<fvVectorMatrix> DivDevTau
        (
            const RhoFieldType& rho,
            volVectorField& U
        ) const;


protected:

    // Protected data

        // Fields

            volSymmTensorField nonlinearStress_;


    // Protected Member Functions

        virtual void correctNonlinearStress(const volTensorField& gradU) = 0;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    // Constructors

        //- Construct from components
        nonlinearEddyViscosity
        (
            const word& modelName,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~nonlinearEddyViscosity()
    {}


    // Member Functions

        //- Return the Reynolds stress tensor [m^2/s^2]
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective surface stress
        virtual tmp<surfaceVectorField> devTau() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonlinearEddyViscosity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
