/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Foam::functionObjects::subtract

Description
    From the first field subtract the remaining fields in the list.

    The operation can be applied to any volume or surface fields generating a
    volume or surface scalar field.

    Example of function object specification:
    \verbatim
    Tdiff
    {
        type            subtract;
        libs            ("libfieldFunctionObjects.so");
        fields          (T Tmean);
        result          Tdiff;
        executeControl  writeTime;
        writeControl    writeTime;
    }
    \endverbatim

See also
    Foam::functionObjects::fieldsExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    subtract.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_subtract_H
#define functionObjects_subtract_H

#include "fieldsExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class subtract Declaration
\*---------------------------------------------------------------------------*/

class subtract
:
    public fieldsExpression
{
    // Private Member Functions

        //- Subtract the list of fields of the specified type
        //  and return the result
        template<class GeoFieldType>
        tmp<GeoFieldType> calcFieldType() const;

        //- Subtract the list of fields and return true if successful
        virtual bool calc();


public:

    friend class fieldsExpression;


    //- Runtime type information
    TypeName("subtract");


    // Constructors

        //- Construct from Time and dictionary
        subtract
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~subtract();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "subtractTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
