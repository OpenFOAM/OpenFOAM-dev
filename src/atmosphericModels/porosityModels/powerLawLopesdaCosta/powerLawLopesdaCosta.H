/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::porosityModels::powerLawLopesdaCosta

Description
    Variant of the power law porosity model with spatially varying
    drag coefficient

    given by:

        \f[
            S = -\rho C_d \Sigma |U|^{(C_1 - 1)} U
        \f]

    where
    \vartable
        \Sigma | Porosity surface area per unit volume
        C_d    | Model linear coefficient
        C_1    | Model exponent coefficient
    \endvartable

    Reference:
    \verbatim
        Costa, J. C. P. L. D. (2007).
        Atmospheric flow over forested and non-forested complex terrain.
    \endverbatim

See also
    Foam::RASModels::kEpsilonLopesdaCosta

SourceFiles
    powerLawLopesdaCosta.C
    powerLawLopesdaCostaTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawLopesdaCosta_H
#define powerLawLopesdaCosta_H

#include "porosityModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                  Class powerLawLopesdaCostaZone Declaration
\*---------------------------------------------------------------------------*/

class powerLawLopesdaCostaZone
{
protected:

    // Protected data

        //- Automatically generated zone name for this porous zone
        const word zoneName_;

        //- Porosity surface area per unit volume zone field
        scalarField Sigma_;


public:

    //- Constructor
    powerLawLopesdaCostaZone
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict
    );

    // Member Functions

        //- Return the porosity surface area per unit volume zone field
        const scalarField& Sigma() const;
};


/*---------------------------------------------------------------------------*\
                      Class powerLawLopesdaCosta Declaration
\*---------------------------------------------------------------------------*/

class powerLawLopesdaCosta
:
    public powerLawLopesdaCostaZone,
    public porosityModel
{
    // Private data

        //- Cd coefficient
        scalar Cd_;

        //- C1 coefficient
        scalar C1_;

        //- Name of density field
        word rhoName_;


    // Private Member Functions

        //- Apply resistance
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Apply resistance
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        powerLawLopesdaCosta(const powerLawLopesdaCosta&);

        //- Disallow default bitwise assignment
        void operator=(const powerLawLopesdaCosta&);


public:

    //- Runtime type information
    TypeName("powerLawLopesdaCosta");

    //- Constructor
    powerLawLopesdaCosta
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~powerLawLopesdaCosta();


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
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;


    // I-O

        //- Write
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "powerLawLopesdaCostaTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
