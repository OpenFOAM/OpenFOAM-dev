/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::fv::option

Description
    Finite volume options abtract base class.  Provides a base set of controls,
    e.g.

        type            scalarExplicitSource    // source type
        active          on;                     // on/off switch
        timeStart       0.0;                    // start time
        duration        1000.0;                 // duration
        selectionMode   cellSet;                // cellSet // points //cellZone
                                                // mapRegion
Note:
    On evaluation, source/sink options are to be added to the equation rhs

SourceFiles
    fvOption.C
    fvOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef fvOption_H
#define fvOption_H

#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"
#include "cellSet.H"
#include "autoPtr.H"
#include "meshToMesh.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class option
{
public:

    // Public data

        //- Enumeration for selection mode types
        enum selectionModeType
        {
            smPoints,
            smCellSet,
            smCellZone,
            smMapRegion,
            smAll
        };

        //- Word list of selection mode type names
        static const NamedEnum<selectionModeType, 5>
            selectionModeTypeNames_;


protected:

    // Protected data

        //- Source name
        const word name_;

        //- Model type
        const word modelType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Top level source dictionary
        dictionary dict_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;

        //- Source active flag
        bool active_;

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- Name of cell set for "cellSet" and "cellZone" selectionMode
        word cellSetName_;

        //- List of points for "points" selectionMode
        List<point> points_;

        //- Set of cells to apply source to
        labelList cells_;

        //- Sum of cell volumes
        scalar V_;

        // Data for smMapRegion only

            //- Mesh to mesh interpolation object
            autoPtr<meshToMesh> meshInterpPtr_;

            //- Name of the neighbour region to map
            word nbrRegionName_;

            //- Master or slave region
            bool master_;


        //- Field names to apply source to - populated by derived models
        wordList fieldNames_;

        //- Applied flag list - corresponds to each fieldNames_ entry
        List<bool> applied_;


    // Protected functions

        //- Flag to bypass the apply flag list checking
        virtual bool alwaysApply() const;

        //- Set the cellSet or points selection
        void setSelection(const dictionary& dict);

        //- Set the cell set based on the user input selection mode
        void setCellSet();


public:

    //- Runtime type information
    TypeName("option");


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            option,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, modelType, dict, mesh)
        );


    // Constructors

        //- Construct from components
        option
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const bool master = false
        );

        //- Return clone
        autoPtr<option> clone() const
        {
            notImplemented("autoPtr<option> clone() const");
            return autoPtr<option>(NULL);
        }

        //- Return pointer to new fvOption object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh database
            const fvMesh& mesh_;
            const word& name_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const word& name
            )
            :
                mesh_(mesh),
                name_(name)
            {}

            autoPtr<option> operator()(Istream& is) const
            {
                //const word name(is);
                const dictionary dict(is);

                return autoPtr<option>
                (
                    option::New(name_, dict, mesh_)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<option> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~option();



    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return dictionary
            inline const dictionary& coeffs() const;

            //- Return const access to the source active flag
            inline bool active() const;

            //- Return const access to the time start
            inline scalar timeStart() const;

            //- Return const access to the duration
            inline scalar duration() const;

            //- Return true if within time limits
            inline bool inTimeLimits(const scalar time) const;

            //- Return const access to the cell selection mode
            inline const selectionModeType& selectionMode() const;

            //- Return const access to the name of cell set for "cellSet"
            //  selectionMode
            inline const word& cellSetName() const;

            //- Return const access to the total cell volume
            inline scalar V() const;

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the mapToMap pointer
            inline const meshToMesh& meshInterp() const;

            //- Return const access to the cell set
            inline const labelList& cells() const;

            //- Set the applied flag to true for field index fieldI
            inline void setApplied(const label fieldI);


        // Edit

            //- Return access to the source active flag
            inline bool& active();

            //- Return access to the time start
            inline scalar& timeStart();

            //- Return access to the duration
            inline scalar& duration();


        // Checks

            //- Is the source active?
            virtual bool isActive();

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;

            //- Check that the source has been applied
            virtual void checkApplied() const;


        // Evaluation

            // Correct

                //- Scalar
                virtual void correct(volScalarField& fld);

                //- Vector
                virtual void correct(volVectorField& fld);

                //- Spherical tensor
                virtual void correct(volSphericalTensorField& fld);

                //- Symmetric tensor
                virtual void correct(volSymmTensorField& fld);

                //- Tensor
                virtual void correct(volTensorField& fld);


            // Add explicit and implicit contributions

                //- Scalar
                virtual void addSup
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Add explicit and implicit contributions to compressible equations

                //- Scalar
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Add explicit and implicit contributions to phase equations

                //- Scalar
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Set values directly

                //- Scalar
                virtual void setValue
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldI
                );

                //- Vector
                virtual void setValue
                (
                    fvMatrix<vector>& eqn,
                    const label fieldI
                );

                //- Spherical tensor
                virtual void setValue
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldI
                );

                //- Symmetric tensor
                virtual void setValue
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldI
                );

                //- Tensor
                virtual void setValue
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldI
                );


            // Flux manipulations

                //- Make the given absolute flux relative
                virtual void makeRelative(surfaceScalarField& phi) const;

                //- Make the given absolute boundary flux relative
                virtual void makeRelative
                (
                    FieldField<fvsPatchField, scalar>& phi
                ) const;

                //- Make the given absolute mass-flux relative
                virtual void makeRelative
                (
                    const surfaceScalarField& rho,
                    surfaceScalarField& phi
                ) const;

                //- Make the given relative flux absolute
                virtual void makeAbsolute(surfaceScalarField& phi) const;

                //- Make the given relative mass-flux absolute
                virtual void makeAbsolute
                (
                    const surfaceScalarField& rho,
                    surfaceScalarField& phi
                ) const;


        // I-O

            //- Write the source header information
            virtual void writeHeader(Ostream&) const;

            //- Write the source footer information
            virtual void writeFooter(Ostream&) const;

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
