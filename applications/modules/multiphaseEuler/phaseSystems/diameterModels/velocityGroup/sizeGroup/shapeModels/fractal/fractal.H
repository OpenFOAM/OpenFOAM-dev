/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

    The effect of sintering, i.e. the coalescence of primary particles
    resulting in a loss of surface area, is taken into account by a separate
    source term in the transport equation for \f$\kappa_i\f$.

Usage
    \table
        Property     | Description             | Required    | Default value
        kappa        | Initial and boundary condition value\\
                                               | yes         |
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

class sinteringModel;

/*---------------------------------------------------------------------------*\
                           Class fractal Declaration
\*---------------------------------------------------------------------------*/

class fractal
:
    public SecondaryPropertyModel<shapeModel>
{
public:

    // Public Data types

        //- Surface growth type enumeration
        enum surfaceGrowthTypes
        {
            sgHardSphere,
            sgParkRogak,
            sgConserved
        };

        //- Surface growth type names
        static const NamedEnum<surfaceGrowthTypes, 3> sgTypeNames_;


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
        volScalarField Su_;

        //- Sintering model
        autoPtr<sinteringModel> sinteringModel_;


    // Private Member Functions

        tmp<volScalarField> dColl() const;


public:

    //- Runtime type information
    TypeName("fractal");


    // Constructors

        //- Construct from dictionary and sizeGroup
        fractal
        (
            const dictionary& dict,
            const sizeGroup& group
        );

        //- Disallow default bitwise copy construction
        fractal(const fractal&) = delete;


    //- Destructor
    virtual ~fractal();


    // Member Functions

        // Access

            //- Return reference to secondary property field
            virtual const volScalarField& fld() const;

            //- Access to secondary property source
            virtual volScalarField& src();

            //- Return representative surface area of the sizeGroup
            virtual const tmp<volScalarField> a() const;

            //- Return representative diameter of the sizeGroup
            virtual const tmp<volScalarField> d() const
            {
                return dColl_;
            }

        // Edit

            //- Correct the collisional diameter
            virtual void correct();

            //- Add drift contribution to secondary property source
            virtual void addDrift
            (
                const volScalarField& Su,
                const sizeGroup& fu,
                const driftModel& model
            );


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
