/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::simplifiedViscousStress

Description
    Simplified linear viscous stress momentum transport model base class

    For incompressible uniform viscosity flow in which the transpose and
    dilatation terms in the linear stress model are formally zero.

    This model can also be useful for complex surface film cases in which
    the transpose and dilatation terms are relatively very small compared to
    the wall surface stress and can safely be ignored avoiding potentially
    large numerical errors around sharp corners.

SourceFiles
    simplifiedViscousStress.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedViscousStress_H
#define simplifiedViscousStress_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class simplifiedViscousStress Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class simplifiedViscousStress
:
    public BasicMomentumTransportModel
{
protected:

    // Protected Member Functions

        //- Return the source term for the momentum equation
        template<class RhoFieldType>
        tmp<fvVectorMatrix> DivDevTau
        (
            const RhoFieldType& rho,
            volVectorField& U
        ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    ClassName("simplified");


    // Constructors

        //- Construct from components
        simplifiedViscousStress
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
    virtual ~simplifiedViscousStress()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read() = 0;

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

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "simplifiedViscousStress.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
