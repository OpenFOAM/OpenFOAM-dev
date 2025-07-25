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
    Foam::fv::fixedValueConstraint

Description
    Constrain the field values within a specified region.

Usage
    For example to set the turbulence properties within a porous region:
    \verbatim
    porosityTurbulence
    {
        type            fixedValueConstraint;

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
#include "fvCellZone.H"
#include "unknownTypeFunction1.H"
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

        //- The cellZone the fvConstraint applies to
        fvCellZone zone_;

        //- Field values
        HashPtrTable<unknownTypeFunction1> fieldValues_;

        //- Fraction of the constraint to apply. Facilitates ramping, or
        //  pulsing, or deactivation after a time. Should take a value between
        //  0 and 1. Defaults to 1 (i.e., apply all of the constraint).
        autoPtr<Function1<scalar>> fraction_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Set value on a field
        template<class Type>
        inline bool constrainType
        (
            fvMatrix<Type>& eqn,
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("fixedValue");


    // Constructors

        //- Construct from components
        fixedValueConstraint
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Add a constraint to an equation
        FOR_ALL_FIELD_TYPES(DEFINE_FV_CONSTRAINT_CONSTRAIN);

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
