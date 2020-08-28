/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::RelativeVelocity

Description
    Function object which generates a lagrangian field containing the velocity
    of each particle relative to the velocity of the surrounding fluid.

    Example usage:
    \verbatim
    relativeVelocity1
    {
        type        relativeVelocity;
    }
    \endverbatim

SourceFiles
    RelativeVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef RelativeVelocity_H
#define RelativeVelocity_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RelativeVelocity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RelativeVelocity
:
    public CloudFunctionObject<CloudType>
{
protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("relativeVelocity");


    // Constructors

        //- Construct from dictionary
        RelativeVelocity
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        RelativeVelocity(const RelativeVelocity<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new RelativeVelocity<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RelativeVelocity();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RelativeVelocity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
