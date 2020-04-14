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
    Foam::LESModels::SpalartAllmarasDDES

Description
    SpalartAllmaras DDES turbulence model for incompressible and compressible
    flows

    Reference:
    \verbatim
        Spalart, P. R., Deck, S., Shur, M. L., Squires, K. D., Strelets, M. K.,
        & Travin, A. (2006).
        A new version of detached-eddy simulation, resistant to ambiguous grid
        densities.
        Theoretical and computational fluid dynamics, 20(3), 181-195.
    \endverbatim

SourceFiles
    SpalartAllmarasDDES.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmarasDDES_H
#define SpalartAllmarasDDES_H

#include "SpalartAllmarasDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                     Class SpalartAllmarasDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class SpalartAllmarasDDES
:
    public SpalartAllmarasDES<BasicMomentumTransportModel>
{
    // Private Member Functions

        tmp<volScalarField::Internal> fd
        (
            const volScalarField::Internal& magGradU
        ) const;

        tmp<volScalarField::Internal> rd
        (
            const volScalarField::Internal& magGradU
        ) const;


protected:

    // Protected Member Functions

        using IOobject::modelName;

        //- Length scale
        virtual tmp<volScalarField::Internal> dTilda
        (
            const volScalarField::Internal& chi,
            const volScalarField::Internal& fv1,
            const volTensorField::Internal& gradU
        ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("SpalartAllmarasDDES");


    // Constructors

        //- Construct from components
        SpalartAllmarasDDES
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
        SpalartAllmarasDDES(const SpalartAllmarasDDES&) = delete;


    //- Destructor
    virtual ~SpalartAllmarasDDES()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SpalartAllmarasDDES&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SpalartAllmarasDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
