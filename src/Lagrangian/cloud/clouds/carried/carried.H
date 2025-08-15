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
    Foam::clouds::carried

Description
    Base class for clouds which are carried by a fluid

SourceFiles
    carried.C

\*---------------------------------------------------------------------------*/

#ifndef carried_H
#define carried_H

#include "cloud.H"
#include "PtrDictionary.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class CarrierField;

namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class carried Declaration
\*---------------------------------------------------------------------------*/

class carried
{
private:

    // Private Data

        //- Carrier fields
        #define DECLARE_CARRIER_FIELDS(Type, nullArg)     \
            mutable PtrDictionary<CarrierField<Type>>     \
                CAT3(carrier, CAPITALIZE(Type), Fields_);
        FOR_ALL_FIELD_TYPES(DECLARE_CARRIER_FIELDS);
        #undef DECLARE_CARRIER_FIELDS

        //- Carrier velocity time derivative
        mutable autoPtr<volVectorField> dUdtcPtr_;


    // Private Member Functions

        //- Access the carrier fields. Only specialisations are defined.
        template<class Type>
        PtrDictionary<CarrierField<Type>>& carrierFields() const;

        //- Clear the carrier fields
        void clearCarrierFields();

        //- Reset the carrier fields
        void resetCarrierFields(const bool predict);

        //- Read the carrier velocity time derivative
        autoPtr<volVectorField> readDUdtc(const cloud& c) const;

        //- Access the carrier velocity time derivative
        const volVectorField& dUdtc() const;


protected:

    // Protected Member Functions

        //- Initialisation hook
        void initialise(const bool predict);

        //- Partition hook
        void partition();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("carried");


    // Public Data

        //- Carrier velocity
        const CarrierField<vector>& Uc;

        //- Carrier velocity curl
        const CarrierField<vector>& curlUc;

        //- Carrier advective derivative
        const CarrierField<vector>& DUDtc;


    // Constructors

        //- Construct from a reference to the cloud
        carried(const cloud& c);


    //- Destructor
    virtual ~carried();


    // Member Functions

        //- Modify a name to disambiguate it as relating to the carrier
        static word carrierName(const word& name);


        // Fields

            //- Add/get a carrier field to/from the cache
            template<class Type>
            const CarrierField<Type>& carrierField(const VolField<Type>&) const;

            //- Add/get a carrier field to/from the cache
            template<class Type, class ... Args>
            const CarrierField<Type>& carrierField(const Args& ... args) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "carriedTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
