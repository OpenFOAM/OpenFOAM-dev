/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::LESModels::kEqn

Description
    One equation eddy-viscosity model

    Eddy viscosity SGS model using a modeled balance equation to simulate the
    behaviour of k.

    Reference:
    \verbatim
        Yoshizawa, A. (1986).
        Statistical theory for compressible turbulent shear flows,
        with the application to subgrid modeling.
        Physics of Fluids (1958-1988), 29(7), 2152-2164.
    \endverbatim

    The default model coefficients are
    \verbatim
        kEqnCoeffs
        {
            Ck                  0.094;
            Ce                  1.048;
        }
    \endverbatim

SourceFiles
    kEqn.C

\*---------------------------------------------------------------------------*/

#ifndef kEqn_H
#define kEqn_H

#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                            Class kEqn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kEqn
:
    public LESeddyViscosity<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Fields

            volScalarField k_;


        // Model constants

            dimensionedScalar Ck_;


    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEqn");


    // Constructors

        //- Constructor from components
        kEqn
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        kEqn(const kEqn&) = delete;


    //- Destructor
    virtual ~kEqn()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                this->nut_ + this->nu()
            );
        }

        //- Correct eddy-Viscosity and related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kEqn&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
