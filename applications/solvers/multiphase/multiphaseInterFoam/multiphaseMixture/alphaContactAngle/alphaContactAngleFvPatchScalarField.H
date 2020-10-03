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
    Foam::alphaContactAngleFvPatchScalarField

Description
    Contact-angle boundary condition for multi-phase interface-capturing
    simulations.  Used in conjunction with multiphaseMixture.

SourceFiles
    alphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaContactAngleFvPatchScalarField_H
#define alphaContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"
#include "multiphaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class alphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class alphaContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    class interfaceThetaProps
    {
        //- Equilibrium contact angle
        scalar theta0_;

        //- Dynamic contact angle velocity scale
        scalar uTheta_;

        //- Limiting advancing contact angle
        scalar thetaA_;

        //- Limiting receding contact angle
        scalar thetaR_;


    public:

        // Constructors
            interfaceThetaProps()
            {}

            interfaceThetaProps(Istream&);


        // Member Functions

            //- Return the equilibrium contact angle theta0
            scalar theta0(bool matched=true) const
            {
                if (matched) return theta0_;
                else return 180.0 - theta0_;
            }

            //- Return the dynamic contact angle velocity scale
            scalar uTheta() const
            {
                return uTheta_;
            }

            //- Return the limiting advancing contact angle
            scalar thetaA(bool matched=true) const
            {
                if (matched) return thetaA_;
                else return 180.0 - thetaA_;
            }

            //- Return the limiting receding contact angle
            scalar thetaR(bool matched=true) const
            {
                if (matched) return thetaR_;
                else return 180.0 - thetaR_;
            }


        // IO functions

            friend Istream& operator>>(Istream&, interfaceThetaProps&);
            friend Ostream& operator<<(Ostream&, const interfaceThetaProps&);
    };

    typedef HashTable
    <
        interfaceThetaProps,
        multiphaseMixture::interfacePair,
        multiphaseMixture::interfacePair::hash
    > thetaPropsTable;


private:

    // Private Data

        thetaPropsTable thetaProps_;


public:

    //- Runtime type information
    TypeName("alphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given alphaContactAngleFvPatchScalarField
        //  onto a new patch
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        alphaContactAngleFvPatchScalarField
        (
            const alphaContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Return the contact angle properties
        const thetaPropsTable& thetaProps() const
        {
            return thetaProps_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
