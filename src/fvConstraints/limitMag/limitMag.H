/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::fv::limitMag

Description
    Limits the magnitude of the specified field to the specified \c max value.

Usage
    Example usage:
    \verbatim
    limitU
    {
        type            limitMag;

        selectionMode   all;

        field           U;

        max             100;
    }
    \endverbatim

SourceFiles
    limitMag.C

\*---------------------------------------------------------------------------*/

#ifndef limitMag_H
#define limitMag_H

#include "fvConstraint.H"
#include "fvCellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class limitMag Declaration
\*---------------------------------------------------------------------------*/

class limitMag
:
    public fvConstraint
{
    // Private data

        //- The set of cells the fvConstraint applies to
        fvCellSet set_;

        //- Field name
        word fieldName_;

        //- Maximum velocity magnitude
        scalar max_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Limit the field
        template<class Type>
        inline bool constrainType
        (
            GeometricField<Type, fvPatchField, volMesh>& psi
        ) const;


public:

    //- Runtime type information
    TypeName("limitMag");


    // Constructors

        //- Construct from components
        limitMag
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        limitMag(const limitMag&) = delete;


    //- Destructor
    virtual ~limitMag()
    {}


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Add a constraint to an equation
        FOR_ALL_FIELD_TYPES(DEFINE_FV_CONSTRAINT_CONSTRAIN_FIELD);

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitMag&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
