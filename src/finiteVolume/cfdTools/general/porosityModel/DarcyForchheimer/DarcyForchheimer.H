/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::porosityModels::DarcyForchheimer

Description
    Darcy-Forchheimer law porosity model, given by:

        \f[
            S = - (\mu d + \frac{\rho |U|}{2} f) U
        \f]

    where
    \vartable
        d        | Darcy coefficient [1/m^2]
        f        | Forchheimer coefficient [1/m]
    \endvartable

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The orientation of the porous region is defined with the same notation as
    a co-ordinate system, but only a Cartesian co-ordinate system is valid.

SourceFiles
    DarcyForchheimer.C
    DarcyForchheimerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef DarcyForchheimer_H
#define DarcyForchheimer_H

#include "porosityModel.H"
#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class DarcyForchheimer Declaration
\*---------------------------------------------------------------------------*/

class DarcyForchheimer
:
    public porosityModel
{
    // Private Data

        //- Darcy coefficient XYZ components (user-supplied) [1/m^2]
        dimensionedVector dXYZ_;

        //- Forchheimer coefficient XYZ components (user-supplied) [1/m]
        dimensionedVector fXYZ_;

        //- Darcy coefficient - converted from dXYZ [1/m^2]
        tensorField D_;

        //- Forchheimer coefficient - converted from fXYZ [1/m]
        tensorField F_;

        //- Name of density field
        word rhoName_;

        //- Name of dynamic viscosity field
        word muName_;

        //- Name of kinematic viscosity field
        word nuName_;


    // Private Member Functions

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;


public:

    //- Runtime type information
    TypeName("DarcyForchheimer");

    // Constructors

        DarcyForchheimer
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& coeffDict,
            const word& cellZoneName
        );

        //- Disallow default bitwise copy construction
        DarcyForchheimer(const DarcyForchheimer&) = delete;


    //- Destructor
    virtual ~DarcyForchheimer();


    // Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData();

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const;

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DarcyForchheimer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DarcyForchheimerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
