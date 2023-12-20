/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::externalCoupledMixedFvPatchField

Description
    This boundary condition provides an interface to an external application.
    Values are transferred as plain text files, where OpenFOAM data is written
    as:

    \verbatim
        # Patch: <patch name>
        <magSf1> <value1> <surfaceNormalGradient1>
        <magSf2> <value2> <surfaceNormalGradient2>
        <magSf3> <value3> <surfaceNormalGradient3>
        ...
        <magSfN> <valueN> <surfaceNormalGradientN>
    \endverbatim

    and received as the constituent pieces of the `mixed' condition, i.e.

    \verbatim
        # Patch: <patch name>
        <value1> <gradient1> <valueFraction1>
        <value2> <gradient2> <valueFraction2>
        <value3> <gradient3> <valueFraction3>
        ...
        <valueN> <gradientN> <valueFractionN>
    \endverbatim

    Data is sent/received as a single file for all patches from the directory

    \verbatim
        $FOAM_CASE/<commsDir>
    \endverbatim

    At start-up, the boundary creates a lock file, i.e..

    \verbatim
        OpenFOAM.lock
    \endverbatim

    ... to signal the external source to wait.  During the boundary condition
    update, boundary values are written to file, e.g.

    \verbatim
        <fileName>.out
    \endverbatim

    The lock file is then removed, instructing the external source to take
    control of the program execution.  When ready, the external program
    should create the return values, e.g. to file

    \verbatim
        <fileName>.in
    \endverbatim

    ... and then re-instate the lock file.  The boundary condition will then
    read the return values, and pass program execution back to OpenFOAM.


Usage
    \table
        Property     | Description             | Required    | Default value
        commsDir     | communications directory   | yes         |
        file         | transfer file name      | yes         |
        waitInterval | interval [s] between file checks | no | 1
        timeOut      | time after which error invoked [s] |no |100*waitInterval
        calcFrequency | calculation frequency  | no          | 1
        initByExternal | external app to initialises values  | yes |
        log          | log program control     | no          | no
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            externalCoupled;
        commsDir        "$FOAM_CASE/comms";
        file            data;
        calcFrequency   1;
        initByExternal  yes;
    }
    \endverbatim

See also
    mixedFvPatchField

SourceFiles
    externalCoupledMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef externalCoupledMixedFvPatchField_H
#define externalCoupledMixedFvPatchField_H

#include "mixedFvPatchFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IFstream;

/*---------------------------------------------------------------------------*\
              Class externalCoupledMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class externalCoupledMixedFvPatchField
:
    public mixedFvPatchField<Type>
{
    // Private Data

        //- Path to communications directory
        fileName commsDir_;

        //- Name of data file
        word fName_;

        //- Interval time between checking for return data [s]
        label waitInterval_;

        //- Time out time [s]
        label timeOut_;

        //- Calculation frequency
        label calcFrequency_;

        //- Flag to indicate values are initialised by external application
        bool initByExternal_;

        //- Log flag
        bool log_;

        //- Master patch flag - controls when to pause/resume execution
        //  Note: only valid when collate option is selected
        bool master_;

        //- Offsets in data file to start reading at correct position
        List<List<label>> offsets_;

        //- Initialised flag
        bool initialised_;

        //- List of coupled patch IDs
        List<label> coupledPatchIndices_;


    // Private Member Functions

        //- Initialise
        void initialise(const fileName& transferFile);

        //- Set the master flag when collate option is selected
        void setMaster(const labelList& patchIDs);

        //- Return the file path to the base communications directory
        fileName baseDir() const;

        //- Write the geometry to the comms dir
        void writeGeometry(OFstream& osPoints, OFstream& osFaces) const;

        //- Return the file path to the lock file
        fileName lockFile() const;

        //- Create lock file
        void createLockFile() const;

        //- Remove lock file
        void removeLockFile() const;

        //- Wait for response from external source
        void startWait() const;

        //- Wait for response from external source
        void wait() const;

        //- Initialise input stream for reading
        void initialiseRead(IFstream& is) const;


protected:

    // Protected Member Functions

        //- Read data from external source
        virtual void readData(const fileName& transferFile);

        //- Write data for external source - calls transferData
        virtual void writeData(const fileName& transferFile) const;

        //- Write header to transfer file
        virtual void writeHeader(OFstream& os) const;


public:

    //- Runtime type information
    TypeName("externalCoupled");

    //- Name of lock file
    static word lockName;

    //- Name of patch key, e.g. '# Patch:' when looking for start of patch data
    static string patchKey;


    // Constructors

        //- Construct from patch, internal field and dictionary
        externalCoupledMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given externalCoupledMixedFvPatchField
        //  onto a new patch
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField&
        ) = delete;

        //- Copy constructor setting internal field reference
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new externalCoupledMixedFvPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~externalCoupledMixedFvPatchField();


    // Member Functions

        //- Return the log flag
        bool log() const
        {
            return log_;
        }

        //- Return the master flag
        bool master() const
        {
            return master_;
        }

        //- Return the master flag
        bool& master()
        {
            return master_;
        }

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Transfer data for external source
        virtual void transferData(OFstream& os) const;

        //- Write the geometry to the comms dir
        void writeGeometry() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "externalCoupledMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
