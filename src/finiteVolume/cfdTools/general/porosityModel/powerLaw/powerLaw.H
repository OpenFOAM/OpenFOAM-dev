/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::powerLaw

Description
    Power law porosity model, given by:

        \f[
            S = - \rho C_0 |U|^{(C_1 - 1)} U
        \f]

    where
    \vartable
        C_0      | model linear coefficient
        C_1      | model exponent coefficient
    \endvartable


SourceFiles
    powerLaw.C
    powerLawTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef powerLaw_H
#define powerLaw_H

#include "porosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    public porosityModel
{
private:

    // Private data

        //- C0 coefficient
        scalar C0_;

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
        powerLaw(const powerLaw&);

        //- Disallow default bitwise assignment
        void operator=(const powerLaw&);


public:

    //- Runtime type information
    TypeName("powerLaw");

    //- Constructor
    powerLaw
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~powerLaw();


    // Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTranformModelData();

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
    #include "powerLawTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
