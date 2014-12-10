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
    Foam::reactingMultiphaseParcelInjectionData

Description
    Container class to provide injection data for reacting multiphase parcels

SourceFiles
    reactingMultiphaseParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef reactingMultiphaseParcelInjectionData_H
#define reactingMultiphaseParcelInjectionData_H

#include "reactingParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class reactingMultiphaseParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const reactingMultiphaseParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    reactingMultiphaseParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
            Class reactingMultiphaseParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class reactingMultiphaseParcelInjectionData
:
    public reactingParcelInjectionData
{
protected:

    // Parcel properties

        //- List of gaseous mass fractions
        scalarList YGas_;

        //- List of liquid mass fractions
        scalarList YLiquid_;

        //- List of solid mass fractions
        scalarList YSolid_;


public:

    //- Runtime type information
    TypeName("reactingMultiphaseParcelInjectionData");

    // Constructors

        //- Null constructor
        reactingMultiphaseParcelInjectionData();

        //- Construct from dictionary
        reactingMultiphaseParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        reactingMultiphaseParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~reactingMultiphaseParcelInjectionData();


    // Access

        //- Return const access to the list of gaseous mass fractions
        inline const scalarList& YGas() const;

        //- Return const access to the list of liquid mass fractions
        inline const scalarList& YLiquid() const;

        //- Return const access to the list of solid mass fractions
        inline const scalarList& YSolid() const;


    // Edit

        //- Return access to the gaseous mass fractions
        inline scalarList& YGas();

        //- Return access to the liquid mass fractions
        inline scalarList& YLiquid();

        //- Return access to the solid mass fractions
        inline scalarList& YSolid();


    // I-O

        // Friend Operators

            friend bool operator==
            (
                const reactingMultiphaseParcelInjectionData& a,
                const reactingMultiphaseParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const reactingMultiphaseParcelInjectionData&, "
                        "const reactingMultiphaseParcelInjectionData&"
                    ")"
                );

                return false;
            }

            friend bool operator!=
            (
                const reactingMultiphaseParcelInjectionData& a,
                const reactingMultiphaseParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const reactingMultiphaseParcelInjectionData&, "
                        "const reactingMultiphaseParcelInjectionData&"
                    ")"
                );

                return false;
            }

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const reactingMultiphaseParcelInjectionData& data
            );

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                reactingMultiphaseParcelInjectionData& data
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactingMultiphaseParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
