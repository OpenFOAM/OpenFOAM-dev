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
    Foam::fv::zeroDimensionalMassSourceBase

Description
    Base class for zero-dimensional mass source models

SourceFiles
    zeroDimensionalMassSource.C

\*---------------------------------------------------------------------------*/

#ifndef zeroDimensionalMassSourceBase_H
#define zeroDimensionalMassSourceBase_H

#include "massSourceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class zeroDimensionalMassSourceBase Declaration
\*---------------------------------------------------------------------------*/

class zeroDimensionalMassSourceBase
:
    public massSourceBase
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

        //- Return the mass flow rate
        virtual scalar massFlowRate() const = 0;


public:

    //- Runtime type information
    TypeName("zeroDimensionalMassSourceBase");


    // Constructors

        //- Construct from explicit source name and mesh
        zeroDimensionalMassSourceBase
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        // Access

            //- Return the cells that the source applies to
            virtual labelUList cells() const;

            //- Return the number of cells that the source applies to
            virtual label nCells() const;

            //- Return the volume of cells that the source applies to
            virtual scalar V() const;


        // Sources

            //- Return the source value
            virtual dimensionedScalar S() const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        //- Update the model
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
