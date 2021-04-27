/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::fv::fixedValueConstraint

Description
    Constrain the field values within a specified region.

Usage
    For example to set the turbulence properties within a porous region:
    \verbatim
    porosityTurbulence
    {
        type            fixedValueConstraint;

        selectionMode   cellZone;
        cellZone        porosity;

        fieldValues
        {
            k           1;
            epsilon     150;
        }
    }
    \endverbatim

SourceFiles
    fixedValueConstraint.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValueConstraint_H
#define fixedValueConstraint_H

#include "fvConstraint.H"
#include "fvCellSet.H"
#include "objectFunction1.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class fixedValueConstraint Declaration
\*---------------------------------------------------------------------------*/

class fixedValueConstraint
:
    public fvConstraint
{
    // Private Member Data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Field values
        HashPtrTable<objectFunction1> fieldValues_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Set value on a field
        template<class Type>
        bool constrainType(fvMatrix<Type>& eqn, const word& fieldName) const;


public:

    //- Runtime type information
    TypeName("fixedValueConstraint");


    // Constructors

        //- Construct from components
        fixedValueConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Add a constraint to an equation
        FOR_ALL_FIELD_TYPES(DEFINE_FV_CONSTRAINT_CONSTRAIN);

        //- Update for mesh changes
        virtual void updateMesh(const mapPolyMesh&);

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
