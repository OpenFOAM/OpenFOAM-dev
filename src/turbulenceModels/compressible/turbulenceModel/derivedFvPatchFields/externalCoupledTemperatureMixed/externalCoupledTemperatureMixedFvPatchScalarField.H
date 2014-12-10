/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::externalCoupledTemperatureMixedFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a temperatue interface to an external
    application.  Values are transferred as plain text files, where OpenFOAM
    data is written as:

        # Patch: <patch name>
        <magSf1> <value1> <qDot1> <htc1>
        <magSf2> <value2> <qDot2> <htc2>
        <magSf3> <value3> <qDot3> <htc2>
        ...
        <magSfN> <valueN> <qDotN> <htcN>

    and received as the constituent pieces of the `mixed' condition, i.e.

        # Patch: <patch name>
        <value1> <gradient1> <valueFracion1>
        <value2> <gradient2> <valueFracion2>
        <value3> <gradient3> <valueFracion3>
        ...
        <valueN> <gradientN> <valueFracionN>

    Data is sent/received as a single file for all patches from the directory

        $FOAM_CASE/<commsDir>

    At start-up, the boundary creates a lock file, i.e..

        OpenFOAM.lock

    ... to signal the external source to wait.  During the boundary condition
    update, boundary values are written to file, e.g.

        <fileName>.out

    The lock file is then removed, instructing the external source to take
    control of the program execution.  When ready, the external program
    should create the return values, e.g. to file

        <fileName>.in

    ... and then re-instate the lock file.  The boundary condition will then
    read the return values, and pass program execution back to OpenFOAM.


    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        commsDir     | communications directory   | yes         |
        fileName     | transfer file name      | yes         |
        waitInterval | interval [s] between file checks | no | 1
        timeOut      | time after which error invoked [s] |no |100*waitInterval
        calcFrequency | calculation frequency  | no          | 1
        log          | log program control     | no          | no
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            externalCoupledTemperature;
        commsDir        "$FOAM_CASE/comms";
        fileName        data;
        calcFrequency   1;
    }
    \endverbatim

SeeAlso
    mixedFvPatchField
    externalCoupledMixedFvPatchField

SourceFiles
    externalCoupledTemperatureMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef externalCoupledTemperatureMixedFvPatchScalarField_H
#define externalCoupledTemperatureMixedFvPatchScalarField_H

#include "externalCoupledMixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IFstream;

/*---------------------------------------------------------------------------*\
     Class externalCoupledTemperatureMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class externalCoupledTemperatureMixedFvPatchScalarField
:
    public externalCoupledMixedFvPatchField<scalar>
{

protected:

    // Protected Member Functions

        //- Write header to transfer file
        virtual void writeHeader(OFstream& os) const;


public:

    //- Runtime type information
    TypeName("externalCoupledTemperature");


    // Constructors

        //- Construct from patch and internal field
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  externalCoupledTemperatureMixedFvPatchScalarField onto a new patch
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const externalCoupledTemperatureMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const externalCoupledTemperatureMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new externalCoupledTemperatureMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        externalCoupledTemperatureMixedFvPatchScalarField
        (
            const externalCoupledTemperatureMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new externalCoupledTemperatureMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~externalCoupledTemperatureMixedFvPatchScalarField();


    // Member functions

        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Transfer data for external source
            virtual void transferData(OFstream& os) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
