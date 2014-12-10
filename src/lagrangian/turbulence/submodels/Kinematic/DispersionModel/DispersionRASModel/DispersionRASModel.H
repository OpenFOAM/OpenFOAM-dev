/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::DispersionRASModel

Description
    Base class for particle dispersion models based on RAS turbulence.

\*---------------------------------------------------------------------------*/

#ifndef DispersionRASModel_H
#define DispersionRASModel_H

#include "DispersionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DispersionRASModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DispersionRASModel
:
    public DispersionModel<CloudType>
{
protected:

    // Protected data

        // Locally cached turbulence fields

            //- Turbulence k
            const volScalarField* kPtr_;

            //- Take ownership of the k field
            bool ownK_;

            //- Turbulence epsilon
            const volScalarField* epsilonPtr_;

            //- Take ownership of the epsilon field
            bool ownEpsilon_;


    // Protected Functions

        //- Return the k field from the turbulence model
        tmp<volScalarField> kModel() const;

        //- Return the epsilon field from the turbulence model
        tmp<volScalarField> epsilonModel() const;


public:

    //- Runtime type information
    TypeName("dispersionRASModel");


    // Constructors

        //- Construct from components
        DispersionRASModel(const dictionary& dict, CloudType& owner);

        //- Construct copy
        DispersionRASModel(DispersionRASModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType> > clone()
        {
            return autoPtr<DispersionModel<CloudType> >
            (
                new DispersionRASModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DispersionRASModel();


    // Member Functions

        //- Cache carrier fields
        virtual void cacheFields(const bool store);


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DispersionRASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
