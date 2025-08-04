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
    Foam::populationBalance::shapeModels::fractal

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
    Foam::populationBalance::shapeModel

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
namespace populationBalance
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

        //- Fractal dimension. Function of representative spherical diameter.
        autoPtr<Function1<scalar>> Df_;

        //- Scaling pre-factor. Function of representative spherical diameter.
        autoPtr<Function1<scalar>> alphaC_;

        //- Groups' ratios of surface area to volume
        PtrList<volScalarField> kappas_;

        //- Groups' collisional diameters
        PtrList<volScalarField> dColls_;

        //- Groups' explicit sources
        PtrList<volScalarField::Internal> Sus_;


    // Private Member Functions

        //- Compute and return the collisional diameter for a group
        tmp<volScalarField> dColl(const label i) const;


public:

    //- Runtime type information
    TypeName("fractal");


    // Constructors

        //- Construct from a dictionary
        fractal
        (
            const dictionary& dict,
            const populationBalanceModel& popBal
        );

        //- Disallow default bitwise copy construction
        fractal(const fractal&) = delete;


    //- Destructor
    virtual ~fractal();


    // Member Functions

        // Access

            //- Return the fractal dimension for a group
            dimensionedScalar Df(const label i) const;

            //- Return the scaling prefactor for a group
            dimensionedScalar alphaC(const label i) const;

            //- Return a reference to a secondary property field
            virtual const volScalarField& fld(const label i) const;

            //- Access a secondary property source
            virtual volScalarField::Internal& src(const label i);

            //- Return the representative surface area of a group
            virtual tmp<volScalarField> a(const label i) const;

            //- Return the representative diameter of a group
            virtual tmp<volScalarField> d(const label i) const;


        // Edit

            //- Solve equations
            virtual void solve();

            //- Correct the shape
            virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fractal&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace shapeModels
} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
