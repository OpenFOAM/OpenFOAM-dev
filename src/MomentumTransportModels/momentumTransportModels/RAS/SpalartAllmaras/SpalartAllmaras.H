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
    Foam::RASModels::SpalartAllmaras

Description
    Spalart-Allmaras one-eqn mixing-length model for incompressible and
    compressible external flows.

    Reference:
    \verbatim
        Spalart, P.R., & Allmaras, S.R. (1994).
        A one-equation turbulence model for aerodynamic flows.
        La Recherche Aerospatiale, 1, 5-21.
    \endverbatim

    The model is implemented without the trip-term and hence the ft2 term is
    not needed.

    It is necessary to limit the Stilda generation term as the model generates
    unphysical results if this term becomes negative which occurs for complex
    flow.  Several approaches have been proposed to limit Stilda but it is not
    clear which is the most appropriate.  Here the limiter proposed by Spalart
    is implemented in which Stilda is clipped at Cs*Omega with the default value
    of Cs = 0.3.

    The default model coefficients are
    \verbatim
        SpalartAllmarasCoeffs
        {
            Cb1         0.1355;
            Cb2         0.622;
            Cw2         0.3;
            Cw3         2.0;
            Cv1         7.1;
            Cs          0.3;
            sigmaNut    0.66666;
            kappa       0.41;
        }
    \endverbatim

SourceFiles
    SpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmaras_H
#define SpalartAllmaras_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class SpalartAllmaras
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
    // Private Member Functions

        using IOobject::modelName;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cs_;


        // Fields

            volScalarField nuTilda_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField::Internal& y_;


    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField::Internal> fv2
        (
            const volScalarField::Internal& chi,
            const volScalarField::Internal& fv1
        ) const;

        tmp<volScalarField::Internal> Stilda
        (
            const volScalarField::Internal& chi,
            const volScalarField::Internal& fv1
        ) const;

        tmp<volScalarField::Internal> fw
        (
            const volScalarField::Internal& Stilda
        ) const;

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- Construct from components
        SpalartAllmaras
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
        SpalartAllmaras(const SpalartAllmaras&) = delete;


    //- Destructor
    virtual ~SpalartAllmaras()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SpalartAllmaras&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmaras.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
