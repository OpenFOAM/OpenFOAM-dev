/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::thermoParcelInjectionData

Description
    Container class to provide injection data for thermodynamic parcels

SourceFiles
    thermoParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef thermoParcelInjectionData_H
#define thermoParcelInjectionData_H

#include "kinematicParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class thermoParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const thermoParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    thermoParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
                 Class thermoParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class thermoParcelInjectionData
:
    public kinematicParcelInjectionData
{
protected:

    // Parcel properties

        //- Temperature [K]
        scalar T_;

        //- Specific heat capacity [J/kg/K]
        scalar Cp_;


public:

    //- Runtime type information
    TypeName("thermoParcelInjectionData");

    // Constructors

        //- Null constructor
        thermoParcelInjectionData();

        //- Construct from dictionary
        thermoParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        thermoParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~thermoParcelInjectionData();


    // Access

        //- Return const access to the temperature
        inline scalar T() const;

        //- Return const access to the specific heat capacity
        inline scalar Cp() const;


    // Edit

        //- Return access to the temperature
        inline scalar& T();

        //- Return access to the specific heat capacity
        inline scalar& Cp();


    // I-O

        // Friend Operators

            friend bool operator==
            (
                const thermoParcelInjectionData& a,
                const thermoParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const thermoParcelInjectionData&, "
                        "const thermoParcelInjectionData&"
                    ")"
                );

                return false;
            }

            friend bool operator!=
            (
                const thermoParcelInjectionData& a,
                const thermoParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const thermoParcelInjectionData&, "
                        "const thermoParcelInjectionData&"
                    ")"
                );

                return false;
            }

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const thermoParcelInjectionData& data
            );

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                thermoParcelInjectionData& data
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
