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
#include "HashPtrTable.H"
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

        //- The name of the carrier phase
        const word carrierPhaseName_;

        //- The name of the corresponding Eulerian phase
        const word phaseName_;

        //- Carrier fields
        #define DECLARE_CARRIER_FIELDS(Type, nullArg)    \
            mutable HashPtrTable<CarrierField<Type>>     \
                CAT3(carrier, CAPITALIZE(Type), Fields_);
        FOR_ALL_FIELD_TYPES(DECLARE_CARRIER_FIELDS);
        #undef DECLARE_CARRIER_FIELDS


    // Private Member Functions

        //- Access the carrier fields. Only specialisations are defined.
        template<class Type>
        HashPtrTable<CarrierField<Type>>& carrierFields() const;

        //- Read the carrier velocity time derivative
        autoPtr<volVectorField> readDUdtc(const cloud& c) const;

        //- Access the carrier velocity time derivative
        const volVectorField& dUdtc() const;


protected:

    // Protected Member Functions

        //- Clear the carrier fields
        void clearCarrierFields();

        //- Reset the carrier fields
        void resetCarrierFields(const bool initial);


public:

    // Public Static Data

        //- Run-time type information
        TypeName("carried");


    // Public Data

        //- Carrier velocity
        const CarrierField<vector>& Uc;

        //- Carrier velocity curl
        const CarrierField<vector>& curlUc;

        //- Carrier velocity time derivative
        private: mutable autoPtr<volVectorField> dUdtcPtr_; public:

        //- Carrier advective derivative
        const CarrierField<vector>& DUDtc;

        //- Corresponding Eulerian phase velocity
        const CarrierField<vector>& UcPhase;


    // Constructors

        //- Construct from a reference to the cloud
        carried(const cloud& c, const dictionary& dict);


    //- Destructor
    virtual ~carried();


    // Member Functions

        //- Return the name of the carrier phase
        const word& carrierPhaseName() const;

        //- Return the name of the corresponding Eulerian phase
        const word& phaseName() const;

        //- Return whether this cloud has a corresponding Eulerian phase
        bool hasPhase() const;

        //- Convert a name to its disambiguated carrier equivalent name. I.e.,
        //  add a 'c' suffix to the member.
        static word nameToCarrierName(const word& name);

        //- Convert a disambiguated carrier name to its equivalent name. I.e.,
        //  remove the 'c' suffix from the member.
        static word carrierNameToName(const word& namec);


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
