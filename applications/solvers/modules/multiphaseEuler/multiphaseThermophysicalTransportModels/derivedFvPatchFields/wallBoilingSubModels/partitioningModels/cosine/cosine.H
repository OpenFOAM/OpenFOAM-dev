/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::wallBoilingModels::partitioningModels::cosine

Description
    Cosine wall heat flux partitioning model.

    Proposed threshold liquid fractions:
      - alphaLiquid1 0.1
      - alphaLiquid0 0.05

    \verbatim
        Tentner, A., Lo, S., & Kozlov, V. (2006).
        Advances in computational fluid dynamics modeling
        of two-phase flow in boiling water reactor fuel assemblies.
        In International Conference of Nuclear Engineering,
        Miami, Florida, USA.
    \endverbatim

SourceFiles
    cosine.C

\*---------------------------------------------------------------------------*/

#ifndef cosine_H
#define cosine_H

#include "partitioningModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace partitioningModels
{

/*---------------------------------------------------------------------------*\
                           Class cosine Declaration
\*---------------------------------------------------------------------------*/

class cosine
:
    public partitioningModel
{
    // Private Data

        //- Lower threshold liquid phase fraction
        scalar alphaLiquid0_;

        //- Upper threshold liquid phase fraction
        scalar alphaLiquid1_;


public:

    //- Runtime type information
    TypeName("cosine");


    // Constructors

        //- Construct from a dictionary
        cosine(const dictionary& dict);

        //- Copy construct
        cosine(const cosine& model);

        //- Construct and return a clone
        virtual autoPtr<partitioningModel> clone() const
        {
            return autoPtr<partitioningModel>(new cosine(*this));
        }


    //- Destructor
    virtual ~cosine();


    // Member Functions

        //- Return the wet fraction
        virtual tmp<scalarField> wetFraction
        (
            const scalarField& alphaLiquid
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace partitioningModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
