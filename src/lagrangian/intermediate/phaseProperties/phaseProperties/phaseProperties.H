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
    Foam::phaseProperties

Description
    Helper class to manage multi-component phase properties

SourceFiles
    phaseProperties.C
    phasePropertiesIO.C

\*---------------------------------------------------------------------------*/

#ifndef phaseProperties_H
#define phaseProperties_H

#include "NamedEnum.H"
#include "Tuple2.H"
#include "PtrList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseProperties Declaration
\*---------------------------------------------------------------------------*/

class phaseProperties
{
public:

    // Public data

        //- Phase type enumeration
        enum phaseType
        {
            GAS,
            LIQUID,
            SOLID,
            UNKNOWN
        };

        //- Corresponding word representations for phase type enumerations
        static const NamedEnum<phaseType, 4> phaseTypeNames_;


private:

   // Private data

        //- Phase type
        phaseType phase_;

        //- State label (s), (l), (g) etc.
        word stateLabel_;

        //- List of component names
        List<word> names_;

        //- List of component mass fractions
        scalarField Y_;

        //- Global ids
        labelList globalIds_;

        //- Map to carrier global id
        labelList globalCarrierIds_;


    // Private Member Functions

        //- Set global ids
        void setGlobalIds(const wordList& globalNames);

        //- Set global carrier ids - attempts to map component names to global
        //  carrier species
        void setGlobalCarrierIds(const wordList& carrierNames);

        //- Check the total mass fraction
        void checkTotalMassFraction() const;

        //- Set the state label
        word phaseToStateLabel(const phaseType pt) const;


public:

    // Constructors

        //- Null constructor
        phaseProperties();

        //- Construct from Istream
        phaseProperties(Istream&);

        //- Construct as copy
        phaseProperties(const phaseProperties&);


    //- Destructor
    ~phaseProperties();


    // Public Member Functions

        //- Initialise the global ids
        void initialiseGlobalIds
        (
            const wordList& gasNames,
            const wordList& liquidNames,
            const wordList& solidNames
        );


        // Access

            //- Return const access to the phase type
            phaseType phase() const;

            //- Return const access to the phase state label
            const word& stateLabel() const;

            //- Return word representation of the phase type
            word phaseTypeName() const;

            //- Return the list of component names
            const List<word>& names() const;

            //- Return const access to a component name
            const word& name(const label cmptI) const;

            //- Return const access to all component mass fractions
            const scalarField& Y() const;

            //- Return non-const access to a component mass fraction
            scalar& Y(const label cmptI);

            //- Return const access to the global ids
            const labelList& globalIds() const;

            //- Return const access to the map to the carrier global ids
            const labelList& globalCarrierIds() const;

            //- Return the global id of a component in the local list by name
            //  Returns -1 if not found
            label globalId(const word& cmptName) const;

            //- Return the id of a component in the local list by name
            //  Returns -1 if not found
            label id(const word& cmptName) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, phaseProperties&);
        friend Ostream& operator<<(Ostream&, const phaseProperties&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
