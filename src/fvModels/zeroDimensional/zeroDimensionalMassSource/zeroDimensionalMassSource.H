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

    Values should be provided for all solved for fields. Warnings will be
    issued if values are not provided for fields for which transport equations
    are solved. Warnings will also be issued if values are provided for fields
    which are not solved for.

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
protected:

    // Protected Data

        //- Cached initial density field
        autoPtr<volScalarField::Internal> rho0Ptr_;


    // Protected Member Functions

        //- Non-virtual read
        void readCoeffs();

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

        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
