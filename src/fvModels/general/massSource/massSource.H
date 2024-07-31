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
    Foam::fv::massSource

Description
    This fvModel applies a mass source to the continuity equation and to all
    field equations. It can be applied to compressible solvers, such as fluid,
    isothermalFluid, compressibleVoF and multiphaseEuler. For incompressible
    solvers, use the volumeSource model instead.

    This model requires a corresponding field source to be specified for all
    solved-for fields.

Usage
    Example usage for a constant mass flow rate applied to a cell set:
    \verbatim
    massSource
    {
        type            massSource;

        select          cellSet;
        cellSet         massSource;

        massFlowRate    1e-4;
    }
    \endverbatim

    Example usage for a pulsing flow rate applied at a point:
    \verbatim
    massSource
    {
        type            massSource;

        select          points;
        points          ((2.75 0.5 0));

        massFlowRate
        {
            type            scale;
            scale           squarePulse;
            start           0.2;
            duration        2;
            value           1e-4;
        }
    }
    \endverbatim

SourceFiles
    massSource.C

See also
    Foam::fv::volumeSource

\*---------------------------------------------------------------------------*/

#ifndef massSource_H
#define massSource_H

#include "massSourceBase.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class massSource Declaration
\*---------------------------------------------------------------------------*/

class massSource
:
    public massSourceBase
{
private:

    // Private Data

        //- The set of cells the source applies to
        autoPtr<fvCellSet> setPtr_;

        //- Mass flow rate
        autoPtr<Function1<scalar>> massFlowRate_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("massSource");


    // Constructors

        //- Construct from explicit source name and mesh
        massSource
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


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
