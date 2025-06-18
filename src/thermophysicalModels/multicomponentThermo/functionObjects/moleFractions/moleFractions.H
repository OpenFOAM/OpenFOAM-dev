/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
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
    Foam::functionObjects::moleFractions

Description
    This function object calculates mole-fraction fields from the mass-fraction
    fields of the multicomponent thermo. The names of the mole-fraction fields
    are obtained from the corresponding mass-fraction fields prepended by "X_".

    Example of function object specification:
    \verbatim
    moleFractions
    {
        type moleFractions;
    }
    \endverbatim

    Optionally, the name of the phase can be specified for multiphase cases.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    moleFractions.C

\*---------------------------------------------------------------------------*/

#ifndef moleFractions_functionObject_H
#define moleFractions_functionObject_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class moleFractions Declaration
\*---------------------------------------------------------------------------*/

class moleFractions
:
    public functionObjects::fvMeshFunctionObject
{
    // Private Data

        //- Optional phase name
        word phaseName_;

        //- Species mole fractions
        PtrList<volScalarField> X_;


public:

    //- Runtime type information
    TypeName("moleFractions");


    // Constructors

        //- Construct from Time and dictionary
        moleFractions(const word& name, const Time& t, const dictionary& dict);

        //- Disallow default bitwise copy construction
        moleFractions(const moleFractions&) = delete;


    //- Destructor
    virtual ~moleFractions();


    // Member Functions

        //- Read the moleFractions data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the mole-fraction fields
        virtual bool execute();

        //- The mole-fraction fields auto-write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const moleFractions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
