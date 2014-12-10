/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::IsotropyModels::NoIsotropy

Description

SourceFiles
    NoIsotropy.C

\*---------------------------------------------------------------------------*/

#ifndef NoIsotropy_H
#define NoIsotropy_H

#include "IsotropyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace IsotropyModels
{

/*---------------------------------------------------------------------------*\
                         Class NoIsotropy Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoIsotropy
:
    public IsotropyModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        NoIsotropy(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoIsotropy(const NoIsotropy<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<IsotropyModel<CloudType> > clone() const
        {
            return autoPtr<IsotropyModel<CloudType> >
            (
                new NoIsotropy<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoIsotropy();


    //- Member Functions

        //- calculate velocities
        virtual void calculate();

        //- Return the model 'active' status
        virtual bool active() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace IsotropyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoIsotropy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
