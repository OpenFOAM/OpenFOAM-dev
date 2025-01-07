/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::LagrangianAverage

Description
    Basic Lagrangian averaging process in which values are averaged in the
    cells and then interpolated assuming a constant value within the cell

SourceFiles
    cell_LagrangianAverage.C
    cell_LagrangianAverages.C

\*---------------------------------------------------------------------------*/

#ifndef cell_LagrangianAverage_H
#define cell_LagrangianAverage_H

#include "LagrangianAverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LagrangianAverages
{

/*---------------------------------------------------------------------------*\
                      Class LagrangianAverage Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cell
:
    public LagrangianAverage<Type>
{
private:

    // Private Classes

        //- Data structure for an average
        struct data
        {
            //- Map from the cell to the cell average
            List<label> cellCellAvg_;

            //- Map from the cell average to the cell
            DynamicList<label> cellAvgCell_;

            //- Number of samples in each cell average
            DynamicList<label> cellAvgCount_;

            //- Weight sums for each cell average
            autoPtr<DynamicList<scalar>> cellAvgWeightSumPtr_;

            //- Value sums for each cell average
            DynamicList<Type> cellAvgSum_;

            //- Construct given a number of cells and a flag to specify whether
            //  or not weight sums need to be stored
            data(const label nCells, const bool hasWeightSum);
        };


    // Private Data

        //- Constant cell weight sums. E.g., cell volume. Might be null.
        const Field<scalar>& weightSum_;

        //- Base data
        data data_;

        //- Flag to indicate whether or not difference data is available
        bool dDataIsValid_;

        //- Difference data for correction later
        data dData_;


    // Private Member Functions

        //- Remove all values from the average data
        static void clear(data& d);

        //- Remove values from the average data
        static void remove
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi,
            data& d
        );

        //- Add values to the average data
        static void add
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi,
            data& d
        );

        //- Interpolate into a sub-field, where possible. Calling code deals
        //  with default values where no elements are available to interpolate.
        virtual void interpolate(LagrangianSubField<Type>& result) const;


public:

    //- Runtime type information
    TypeName("cell");


    // Constructors

        //- Construct with a name, for a mesh and with given dimensions
        cell
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dimensionSet& dimensions,
            const Field<scalar>& weightSum
        );


    //- Destructor
    virtual ~cell();


    // Member Functions

        //- Remove weighted values from the average
        virtual void remove
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi
        );

        //- Add weighted values to the average
        virtual void add
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi,
            const bool cache
        );

        //- Correct weighted values in the average
        virtual void correct
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi,
            const bool cache
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LagrangianAverages
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cell_LagrangianAverage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
