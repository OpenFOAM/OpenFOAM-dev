/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

    This model requires a corresponding field source to be specified for all
    solved-for fields.

    This model will write out additional fields, zeroDimensionalMassSource:m
    and zeroDimensionalMassSource:factor. The zeroDimensionalMassSource:m field
    is the total accumulated mass; the sum of the starting mass, plus all mass
    added by models of this type. This may differ from the actual current mass
    (= rho*V) within the zero-dimensional cell if pressure or density
    constraints are being used. The zeroDimensionalMassSource:factor field is
    the ratio between the current mass or volume and the total accumulated mass
    or volume.

Usage
    Example usage:
    \verbatim
    zeroDimensionalMassSource
    {
        type            zeroDimensionalMassSource;

        massFlowRate    1e-4;
    }
    \endverbatim

SourceFiles
    zeroDimensionalMassSource.C

\*---------------------------------------------------------------------------*/

#ifndef zeroDimensionalMassSource_H
#define zeroDimensionalMassSource_H

#include "zeroDimensionalMassSourceBase.H"
#include "Function1.H"

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
    public zeroDimensionalMassSourceBase
{
private:

    // Private Data

        //- Zero-dimensional mass flow rate
        autoPtr<Function1<scalar>> zeroDimensionalMassFlowRate_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

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

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
