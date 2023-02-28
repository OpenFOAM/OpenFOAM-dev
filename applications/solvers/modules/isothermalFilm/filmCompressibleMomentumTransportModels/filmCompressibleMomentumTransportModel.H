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
    Foam::filmCompressibleMomentumTransportModel

Description
    Abstract base class for film compressible momentum transport models.

SourceFiles
    filmCompressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef filmCompressibleMomentumTransportModel_H
#define filmCompressibleMomentumTransportModel_H

#include "compressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class filmCompressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class filmCompressibleMomentumTransportModel
:
    public compressibleMomentumTransportModel
{

public:

    typedef volScalarField alphaField;
    typedef volScalarField rhoField;


protected:

    // Protected data

        const alphaField& alpha_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            filmCompressibleMomentumTransportModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const viscosity& viscosity
            ),
            (alpha, rho, U, alphaRhoPhi, phi, viscosity)
        );


    // Constructors

        //- Construct
        filmCompressibleMomentumTransportModel
        (
            const word& type,
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<filmCompressibleMomentumTransportModel> New
        (
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~filmCompressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devTau() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const = 0;
};


namespace filmCompressible
{
    typedef filmCompressibleMomentumTransportModel momentumTransportModel;

    template<class FilmCompressibleMomentumTransportModel>
    autoPtr<FilmCompressibleMomentumTransportModel> New
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "filmCompressibleMomentumTransportModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
