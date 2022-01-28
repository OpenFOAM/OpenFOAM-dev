/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    simulations.

SourceFiles
    alphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphaContactAngleFvPatchScalarField_H
#define alphaContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class alphaContactAngleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphaContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
public:

    class contactAngleProperties
    {
        // Private Data

            //- Equilibrium contact angle
            scalar theta0_;

            //- Is this contact angle dynamic?
            bool dynamic_;

            //- Dynamic contact angle velocity scale
            scalar uTheta_;

            //- Limiting advancing contact angle
            scalar thetaA_;

            //- Limiting receding contact angle
            scalar thetaR_;


    public:

        // Constructors

            //- Construct null
            contactAngleProperties();

            //- Construct non-dynamic properties from components
            contactAngleProperties(const scalar theta0);

            //- Construct dynamic properties from components
            contactAngleProperties
            (
                const scalar theta0,
                const scalar uTheta,
                const scalar thetaA,
                const scalar thetaR
            );

            //- Construct from a dictionary
            contactAngleProperties(const dictionary& dict);


        // Member Functions

            // Access

                //- Return the equilibrium contact angle theta0
                inline scalar theta0() const
                {
                    return theta0_;
                }

                //- Return if this contact angle is dynamic
                inline bool dynamic() const
                {
                    return dynamic_;
                }

                //- Return the dynamic contact angle velocity scale
                inline scalar uTheta() const
                {
                    return uTheta_;
                }

                //- Return the limiting advancing contact angle
                inline scalar thetaA() const
                {
                    return thetaA_;
                }

                //- Return the limiting receding contact angle
                inline scalar thetaR() const
                {
                    return thetaR_;
                }


            //- Return the contact angle properties for the reverse of this
            //  interface (e.g., convert air_water to water_air).
            contactAngleProperties reversed() const;

            //- Write to stream
            void write(Ostream& os) const;


        // Member operators

            //- Check equality between two sets of contact angle properties
            bool operator==(const contactAngleProperties& thetaProps) const;

            //- Check inequality between two sets of contact angle properties
            bool operator!=(const contactAngleProperties& thetaProps) const;
    };


private:

    // Private Data

        //- Table of contact angle properties
        HashTable<contactAngleProperties> thetaProps_;


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
        const HashTable<contactAngleProperties>& thetaProps() const
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
