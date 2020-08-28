/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::laminarModels::Giesekus

Description
    Giesekus model for viscoelasticity using the upper-convected time
    derivative of the stress tensor with support for multiple modes.

    Reference:
    \verbatim
        Giesekus, H., 1982.
        A simple constitutive equation for polymer fluids based on the
        concept of deformation-dependent tensional mobility.
        J. Non-Newton. Fluid. 11, 69–109.
    \endverbatim

See also
    Foam::laminarModels::Maxwell

SourceFiles
    Giesekus.C

\*---------------------------------------------------------------------------*/

#ifndef Giesekus_H
#define Giesekus_H

#include "Maxwell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

/*---------------------------------------------------------------------------*\
                          Class Giesekus Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class Giesekus
:
    public Maxwell<BasicMomentumTransportModel>
{
protected:

    // Protected data

        // Mode coefficients

            PtrList<dimensionedScalar> alphaGs_;


    // Protected Member Functions

        virtual tmp<fvSymmTensorMatrix> sigmaSource
        (
            const label modei,
            volSymmTensorField& sigma
        ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("Giesekus");


    // Constructors

        //- Construct from components
        Giesekus
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
        Giesekus(const Giesekus&) = delete;



    //- Destructor
    virtual ~Giesekus()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Giesekus&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Giesekus.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
