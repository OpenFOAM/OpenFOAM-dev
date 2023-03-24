/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fv::zeroDimensionalMassSource

Description
    This fvModel applies a mass source to the continuity equation and to all
    field equations, in a zero-dimensional case. Correction is made to account
    for the mass that exits the domain due to expansion in space, so that the
    model correctly applies a total mass flow rate.

    This model will write out additional fields, zeroDimensionalMassSource:m
    and zeroDimensionalMassSource:factor. The zeroDimensionalMassSource:m field
    is the total accumulated mass; the sum of the starting mass, plus all mass
    added by models of this type. This may differ from the actual current mass
    (= rho*V) within the zero-dimensional cell if pressure or density
    constraints are being used. The zeroDimensionalMassSource:factor field is
    the ratio between the current mass or volume and the total accumulated mass
    or volume.

    If the mass flow rate is positive then user-supplied fixed property values
    are introduced to the field equations. If the mass flow rate is negative
    then properties are removed at their current value.

Usage
    Example usage:
    \verbatim
    zeroDimensionalMassSource
    {
        type            zeroDimensionalMassSource;

        massFlowRate    1e-4;

        fieldValues
        {
            U               (10 0 0);
            T               300;
            k               0.375;
            epsilon         14.855;
        }
    }
    \endverbatim

    If the mass flow rate is positive then values should be provided for all
    solved for fields. Warnings will be issued if values are not provided for
    fields for which transport equations are solved. Warnings will also be
    issued if values are provided for fields which are not solved for.

SourceFiles
    zeroDimensionalMassSource.C

\*---------------------------------------------------------------------------*/

#ifndef zeroDimensionalMassSource_H
#define zeroDimensionalMassSource_H

#include "massSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class zeroDimensionalMassSource Declaration
\*---------------------------------------------------------------------------*/

class zeroDimensionalMassSource
:
    public massSource
{
private:

    // Private Data

        //- Reference to the zero-dimensional mass
        volScalarField& m0D_;


    // Private Member Functions

        //- Calculate and return the zero-dimensional mass
        tmp<volScalarField> calcM0D() const;

        //- Initialise and return a reference to the zero-dimensional mass
        volScalarField& initM0D() const;

        //- Get a reference to the mass
        const volScalarField& m() const;


protected:

    // Protected Member Functions

        //- Return the mass flow rate
        virtual scalar massFlowRate() const;


public:

    //- Runtime type information
    TypeName("zeroDimensionalMassSource");


    // Constructors

        //- Construct from explicit source name and mesh
        zeroDimensionalMassSource
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        //- Update the model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
