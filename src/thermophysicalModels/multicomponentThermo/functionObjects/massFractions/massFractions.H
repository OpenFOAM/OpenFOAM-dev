/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::functionObjects::massFractions

Description
    This function object calculates mass-fraction fields from mole-fraction or
    moles fields present on disk. This is intended to be used for
    initialisation where mole-fractions are known. If any mass fraction fields
    are found (other than Ydefault) then an error will be generated and the
    fields will not be overwritten. The names of the mole-fraction fields are
    obtained from the corresponding mass-fraction fields prepended by "X_", and
    the moles fields are prepended by "n_". Either mole-fraction fields or
    moles fields should be present, not both.

    Example of function object specification:
    \verbatim
    massFractions
    {
        type massFractions;
    }
    \endverbatim

    Optionally, the name of the phase can be specified for multiphase cases.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    massFractions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_massFractions_H
#define functionObjects_massFractions_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class massFractions Declaration
\*---------------------------------------------------------------------------*/

class massFractions
:
    public functionObjects::fvMeshFunctionObject
{
    // Private Data

        //- Optional phase name
        word phaseName_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


public:

    //- Runtime type information
    TypeName("massFractions");


    // Constructors

        //- Construct from Time and dictionary
        massFractions(const word& name, const Time& t, const dictionary& dict);

        //- Disallow default bitwise copy construction
        massFractions(const massFractions&) = delete;


    //- Destructor
    virtual ~massFractions();


    // Member Functions

        //- Read the massFractions data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the mass-fraction fields
        virtual bool execute();

        //- The mass-fraction fields auto-write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const massFractions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
