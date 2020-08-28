/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::distributionModels::multiNormal

Description
    A multiNormal distribution model

    \verbatim
        model = sum_i strength_i * exp(-0.5*((x - expectation_i)/variance_i)^2 )
    \endverbatim


SourceFiles
    multiNormal.C

\*---------------------------------------------------------------------------*/

#ifndef multiNormal_H
#define multiNormal_H

#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                         Class multiNormal Declaration
\*---------------------------------------------------------------------------*/

class multiNormal
:
    public distributionModel
{
    // Private Data

        //- Distribution minimum
        scalar minValue_;

        //- Distribution maximum
        scalar maxValue_;

        //- Distribution range
        scalar range_;


        // Model coefficients

            List<scalar> expectation_;
            List<scalar> variance_;
            List<scalar> strength_;


public:

    //- Runtime type information
    TypeName("multiNormal");


    // Constructors

        //- Construct from components
        multiNormal(const dictionary& dict, Random& rndGen);

        //- Construct copy
        multiNormal(const multiNormal& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new multiNormal(*this));
        }


    //- Destructor
    virtual ~multiNormal();


    // Member Functions

        //- Sample the distributionModel
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;

        //- Return the mean value
        virtual scalar meanValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
