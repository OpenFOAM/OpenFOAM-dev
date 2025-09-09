/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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
    Foam::diameterModels::shapeModels::fractal

Description
    Class for modelling the shape of particle aggregates using the concept of
    fractal geometry. Returns a collisional diameter

    \f[
        d_{c_i} =
            \frac{6}{\kappa_i}
            \left(
                \frac{v_i \kappa_i^3}{36 \pi \alpha_c}
            \right)^{1/D_{f_i}}\,,
    \f]

    computed from a constant fractal dimension \f$D_{f_i}\f$ and a
    field-dependent surface area to volume ratio \f$\kappa_i\f$, assuming that
    the primary particles in an aggregate have the same size.

Usage
    \table
        Property     | Description             | Required    | Default value
        Df           | Fractal dimension       | yes         |
        alphaC       | Scaling prefactor       | yes         |
    \endtable

See also
    Foam::diameterModels::shapeModel

SourceFiles
    fractal.C

\*---------------------------------------------------------------------------*/

#ifndef fractal_H
#define fractal_H

#include "SecondaryPropertyModel.H"
#include "shapeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
namespace shapeModels
{

/*---------------------------------------------------------------------------*\
                           Class fractal Declaration
\*---------------------------------------------------------------------------*/

class fractal
:
    public SecondaryPropertyModel<shapeModel>
{
private:

    // Private Data

        //- Ratio of surface area to volume
        volScalarField kappa_;

        //- Fractal dimension
        dimensionedScalar Df_;

        //- Scaling prefactor
        dimensionedScalar alphaC_;

        //- Collisional diameter
        volScalarField dColl_;

        //- Explicit source
        volScalarField::Internal Su_;


    // Private Member Functions

        //- Compute and return the collisional diameter
        tmp<volScalarField> dColl() const;


public:

    //- Runtime type information
    TypeName("fractal");


    // Constructors

        //- Construct from dictionaries and sizeGroup
        fractal
        (
            const dictionary& dict,
            const sizeGroup& group,
            const dictionary& groupDict
        );

        //- Disallow default bitwise copy construction
        fractal(const fractal&) = delete;


    //- Destructor
    virtual ~fractal();


    // Member Functions

        // Access

            //- Return the fractal dimension
            inline const dimensionedScalar& Df() const
            {
                return Df_;
            }

            //- Return the scaling prefactor
            inline const dimensionedScalar& alphaC() const
            {
                return alphaC_;
            }

            //- Return reference to secondary property field
            virtual const volScalarField& fld() const;

            //- Access to secondary property source
            virtual volScalarField::Internal& src();

            //- Return representative surface area of the sizeGroup
            virtual const tmp<volScalarField> a() const;

            //- Return representative diameter of the sizeGroup
            virtual const tmp<volScalarField> d() const;


        // Edit

            //- Correct the collisional diameter
            virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fractal&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace shapeModels
} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
