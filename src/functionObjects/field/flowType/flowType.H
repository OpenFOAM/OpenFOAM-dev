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
    Foam::functionObjects::flowType

Description
    Calculates and writes a field in the range -1 to 1, in which:
      -  1 = rotational flow
      -  0 = parallel shear (50% rotational, 50% shear)
      - -1 = pure shear deformation

    The field is calculated as:
    \verbatim
    flowType = -(devGradU && devGradU.T())/magSqr(devGradU)
    \endverbatim

    where \c devGradU = \c dev(fvc::grad(U)).  The \c flowType parameter is
    designed for visualisation of vortices, by extracting an iso-surface of
    these fields at a specified value.  The parameter is similar to \c Q and
    \c Lambda2, except that \c Q and \c Lambda2 represent the rotational speed,
    or "strength", of vortices whereas \c flowType is normalised to a scale
    between -1 and 1, making it independent of vortex strength.

    The normalisation of \c flowType provides the advantage that it is much
    easier to select a value of iso-surface to visualise vortices effectively.
    Generally a value of 0.2 (or between 0.05 - 0.2) is recommended.  With \c Q
    and \c Lambda2 the choice of an effective iso-surface value depends on the
    strength of vortices for the specific CFD problem.

    The \c flowType function object writes an additional \c magSqrDevGradU field
    corresponding to the denominator in the \c flowType field.  An iso-surface
    can then be coloured by the \c magSqrDevGradU field to indicate the strength
    of the vortices as they move with the flow.

See also
    Foam::functionObjects::Q
    Foam::functionObjects::Lambda2

SourceFiles
    flowType.C

\*---------------------------------------------------------------------------*/

#ifndef flowType_functionObject_H
#define flowType_functionObject_H

#include "fvMeshFunctionObject.H"
#include "writeLocalObjects.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class flowType Declaration
\*---------------------------------------------------------------------------*/

class flowType
:
    public fvMeshFunctionObject,
    public writeLocalObjects
{
    // Private Member data

        //- Name of the velocity field, default "U"
        const word fieldName_;

        //- Name of the flowType field
        const word resultName_;


public:

    //- Runtime type information
    TypeName("flowType");


    // Constructors

        //- Construct from Time and dictionary
        flowType
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~flowType();

    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate the flowType fields
        virtual bool execute();

        //- Write the flowType fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
