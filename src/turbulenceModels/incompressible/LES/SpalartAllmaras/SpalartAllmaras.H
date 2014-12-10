/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::incompressible::LESModels::SpalartAllmaras

Group
    grpIcoDESTurbulence

Description
    SpalartAllmaras DES (SA + LES) turbulence model for incompressible flows

    Extended according to
    \verbatim
        "An Unstructured Grid Generation and Adaptive Solution Technique
        for High Reynolds Number Compressible Flows"
        G.A. Ashford,
        Ph.D. thesis, University of Michigan, 1996.
    \endverbatim
    by using the optional flag \c ashfordCorrection

SourceFiles
    SpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmaras_H
#define SpalartAllmaras_H

#include "DESModel.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public DESModel
{
    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields();

        // Disallow default bitwise copy construct and assignment
        SpalartAllmaras(const SpalartAllmaras&);
        SpalartAllmaras& operator=(const SpalartAllmaras&);


protected:

    // Protected data

        dimensionedScalar sigmaNut_;
        dimensionedScalar kappa_;


        // Model constants

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;
            dimensionedScalar CDES_;
            dimensionedScalar ck_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;


        //- Optional flag to activate the Ashford correction
        Switch ashfordCorrection_;


        // Fields

            wallDist y_;
            volScalarField nuTilda_;
            volScalarField nuSgs_;


    // Protected Member Functions

        virtual tmp<volScalarField> fv1() const;
        virtual tmp<volScalarField> fv2() const;
        virtual tmp<volScalarField> fv3() const;
        virtual tmp<volScalarField> S(const volTensorField& gradU) const;

        virtual tmp<volScalarField> STilda
        (
            const volScalarField& S,
            const volScalarField& dTilda
        ) const;

        virtual tmp<volScalarField> r
        (
            const volScalarField& visc,
            const volScalarField& S,
            const volScalarField& dTilda
        ) const;

        virtual tmp<volScalarField> fw
        (
            const volScalarField& S,
            const volScalarField& dTilda
        ) const;

        //- Length scale
        virtual tmp<volScalarField> dTilda(const volScalarField& S) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- Construct from components
        SpalartAllmaras
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~SpalartAllmaras()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        tmp<volScalarField> nuTilda() const
        {
            return nuTilda_;
        }

        //- Return SGS viscosity
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Correct nuTilda and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
