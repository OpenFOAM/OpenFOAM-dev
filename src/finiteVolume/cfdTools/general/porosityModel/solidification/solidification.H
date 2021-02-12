/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::porosityModels::solidification

Description
    Simple solidification porosity model

    This is a simple approximation to solidification where the solid phase
    is represented as a porous blockage with the drag-coefficient evaluated from

        \f[
            S = - \alpha \rho D(T) U
        \f]

    where
    \vartable
        \alpha  | Optional phase-fraction of solidifying phase
        D(T)    | User-defined drag-coefficient as function of temperature
    \endvartable

    Note that the latent heat of solidification is not included and the
    temperature is unchanged by the modelled change of phase.

    Example of the solidification model specification:
    \verbatim
        type            solidification;

        solidificationCoeffs
        {
            // Solidify between 330K and 330.5K
            D table
            (
                (330.0     10000) // Solid below 330K
                (330.5     0)     // Liquid above 330.5K
            );

            // Optional phase-fraction of solidifying phase
            alpha alpha.liquid;

            // Solidification porosity is isotropic
            // use the global coordinate system
            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1      (1 0 0);
                    e2      (0 1 0);
                }
            }
        }
    \endverbatim

SourceFiles
    solidification.C
    solidificationTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solidification_H
#define solidification_H

#include "porosityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                       Class solidification Declaration
\*---------------------------------------------------------------------------*/

class solidification
:
    public porosityModel
{
    // Private Data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Name of optional phase-fraction field, default = "none"
        word alphaName_;

        //- Name of density field, default = "rho"
        word rhoName_;

        //- User-defined drag-coefficient as function of temperature
        autoPtr<Function1<scalar>> D_;


    // Private Member Functions

        //- Apply resistance
        template<class AlphaFieldType, class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            const scalarField& V,
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            const volVectorField& U
        ) const;

        //- Apply resistance
        template<class AlphaFieldType, class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            const volVectorField& U
        ) const;

        //- Apply resistance
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            const scalarField& V,
            const RhoFieldType& rho,
            const volVectorField& U
        ) const;

        //- Apply resistance
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const volVectorField& U
        ) const;


public:

    //- Runtime type information
    TypeName("solidification");

    // Constructors

        solidification
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName
        );

        //- Disallow default bitwise copy construction
        solidification(const solidification&) = delete;


    //- Destructor
    virtual ~solidification();


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


    // I-O

        //- Write
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidification&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidificationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
