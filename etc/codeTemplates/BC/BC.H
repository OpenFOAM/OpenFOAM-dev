/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::CONSTRUCT

Description
    This boundary condition provides a NAME condition,
    calculated as:

        \f[
            Q = Q_{0} + Q_{p} + s*Q_{t}
        \f]

    where
    \vartable
        s      | single scalar value [units]
        Q_{0}  | single TYPE value [units]
        Q_{p}  | TYPE field across patch [units]
        Q_{t}  | TYPE function of time [units]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        scalarData | single scalar value          | yes    |
        data       | single TYPE value          | yes    |
        fieldData  | TYPE field across patch    | yes    |
        timeVsData | TYPE function of time      | yes    |
        wordData   | word, eg name of data object | no     | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       NAME;
        scalarData -1;
        data       ONE;
        fieldData  uniform THREE;
        timeVsData table (
                             (0 ZERO)
                             (1 TWO)
                          );
        wordName   anotherName;
        value      uniform FOUR; // optional initial value
    }
    \endverbatim

SourceFiles
    CONSTRUCT.C

\*---------------------------------------------------------------------------*/

#ifndef CONSTRUCT_H
#define CONSTRUCT_H

#include "BASEFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class CONSTRUCT Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CONSTRUCT
:
    public PARENT
{
    // Private data

        //- Single valued scalar quantity, e.g. a coefficient
        scalar scalarData_;

        //- Single valued Type quantity, e.g. reference pressure pRefValue_
        //  Other options include vector, tensor
        TYPE data_;

        //- Field of Types, typically defined across patch faces
        //  e.g. total pressure p0_.  Other options include vectorField
        FIELD fieldData_;

        //- Type specified as a function of time for time-varying BCs
        autoPtr<Function1<TYPE>> timeVsData_;

        //- Word entry, e.g. pName_ for name of the pressure field on database
        word wordData_;

        //- Label, e.g. patch index, current time index
        label labelData_;

        //- Boolean for true/false, e.g. specify if flow rate is volumetric_
        bool boolData_;


    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("NAME");


    // Constructors

        //- Construct from patch and internal field
        CONSTRUCT
        (
            const fvPatch&,
            const DimensionedField<TYPE, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CONSTRUCT
        (
            const fvPatch&,
            const DimensionedField<TYPE, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given BASETypeFvPatchField
        //  onto a new patch
        CONSTRUCT
        (
            const CLASS&,
            const fvPatch&,
            const DimensionedField<TYPE, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        CONSTRUCT
        (
            const CLASS&
        );

        //- Construct and return a clone
        virtual tmp<FVPATCHF> clone() const
        {
            return tmp<FVPATCHF>
            (
                new CLASS(*this)
            );
        }

        //- Construct as copy setting internal field reference
        CONSTRUCT
        (
            const CLASS&,
            const DimensionedField<TYPE, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<FVPATCHF> clone
        (
            const DimensionedField<TYPE, volMesh>& iF
        ) const
        {
            return tmp<FVPATCHF>
            (
                new CLASS
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const FVPATCHF&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CONSTRUCT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
