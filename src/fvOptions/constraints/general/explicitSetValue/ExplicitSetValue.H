/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::fv::explicitSetValue

Description
    Set values field values explicity.

    Sources described by:

    \verbatim
    <Type>ExplicitSetValueCoeffs
    {
        injectionRate
        {
            k           30.7;
            epsilon     1.5;
        }
    }
    \endverbatim

SeeAlso
    Foam::fvOption

SourceFiles
    explicitSetValue.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitSetValue_H
#define ExplicitSetValue_H

#include "fvOption.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class explicitSetValue Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ExplicitSetValue
:
    public option
{

protected:

    // Protected data

        //- Source value per field
        List<Type> injectionRate_;


    // Protected functions

        //- Set the local field data
        void setFieldData(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("explicitSetValue");


    // Constructors

        //- Construct from components
        ExplicitSetValue
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluation

            //- Set value on field
            virtual void setValue(fvMatrix<Type>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ExplicitSetValue.C"
#   include "ExplicitSetValueIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
