/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::PackingModels::NoPacking

Description

SourceFiles
    NoPacking.C

\*---------------------------------------------------------------------------*/

#ifndef NoPacking_H
#define NoPacking_H

#include "PackingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PackingModels
{

/*---------------------------------------------------------------------------*\
                          Class NoPacking Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPacking
:
    public PackingModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoPacking(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoPacking(const NoPacking<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<PackingModel<CloudType>> clone() const
        {
            return autoPtr<PackingModel<CloudType>>
            (
                new NoPacking<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoPacking();


    // Member Functions

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            typename CloudType::parcelType& p,
            const scalar deltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PackingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoPacking.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
