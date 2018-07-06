/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::CloudSubModelBase

Description
    Base class for cloud sub-models

SourceFiles
    CloudSubModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef CloudSubModelBase_H
#define CloudSubModelBase_H

#include "subModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CloudSubModelBase Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudSubModelBase
:
    public subModelBase
{
protected:

    // Protected Data

        //- Reference to the cloud
        CloudType& owner_;


public:

    // Constructors

        //- Construct null from owner cloud
        CloudSubModelBase(CloudType& owner);

        //- Construct from owner cloud without name
        CloudSubModelBase
        (
            CloudType& owner,
            const dictionary& dict,
            const word& baseName,
            const word& modelType,
            const word& dictExt = "Coeffs"
        );

        //- Construct from owner cloud with name
        CloudSubModelBase
        (
            const word& modelName,
            CloudType& owner,
            const dictionary& dict,
            const word& baseName,
            const word& modelType
        );

        //- Construct as copy
        CloudSubModelBase(const CloudSubModelBase<CloudType>& smb);


    //- Destructor
    virtual ~CloudSubModelBase();

    //- Type of cloud this model was instantiated for
    typedef CloudType cloudType;


    // Member Functions

        // Access

            //- Return const access to the owner cloud
            const CloudType& owner() const;

            //- Flag to indicate when to write a property
            virtual bool writeTime() const;


        // Edit

            //- Return non-const access to the owner cloud for manipulation
            CloudType& owner();


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudSubModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
