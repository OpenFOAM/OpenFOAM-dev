/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "tecplotWriter.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::tecplotWriter::tecplotWriter(const Time& runTime)
:
    runTime_(runTime)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tecplotWriter::writeInit
(
    const word& name,
    const string& varNames,
    const fileName& fName,
    INTEGER4 tecplotFileType
) const
{
Pout<< endl
    << endl
    << "Name:" << name
    << " varNames:" << varNames
    << " to file:" << fName
    << " of type:" << tecplotFileType
    << endl;

    INTEGER4 IsDouble = 0;  // float
    INTEGER4 Debug = 0;     // nodebug
    if
    (
        !TECINI112
        (
            const_cast<char*>(name.c_str()),       /* Data Set Title       */
            const_cast<char*>(varNames.c_str()),   /* Variable List        */
            const_cast<char*>(fName.c_str()),      /* File Name            */
            const_cast<char*>(runTime_.path().c_str()), /* Scratch Directory */
            &tecplotFileType,
            &Debug,
            &IsDouble
        )
    )
    {
//        FatalErrorInFunction
//            << "Error in TECINI112." << exit(FatalError);
    }
}


void Foam::tecplotWriter::writePolyhedralZone
(
    const word& zoneName,
    INTEGER4 strandID,
    const fvMesh& mesh,
    const List<INTEGER4>& varLocArray,
    INTEGER4 nFaceNodes
) const
{
    /* Call TECZNE112 */
    INTEGER4  NumNodes   = mesh.nPoints();         /* number of unique nodes */
    INTEGER4  NumElems   = mesh.nCells();         /* number of elements */
    INTEGER4  NumFaces   = mesh.nFaces();         /* number of unique faces */

    INTEGER4  ICellMax   = 0;         /* Not Used, set to zero */
    INTEGER4  JCellMax   = 0;         /* Not Used, set to zero */
    INTEGER4  KCellMax   = 0;         /* Not Used, set to zero */

    double    SolTime    = runTime_.value();     /* solution time   */
    INTEGER4  ParentZone = 0;         /* no parent zone  */

    INTEGER4  IsBlock    = 1;         /* block format  */

    INTEGER4  NFConns    = 0;         /* not used for FEPolyhedron
                                       * zones
                                       */
    INTEGER4  FNMode     = 0;         /* not used for FEPolyhedron
                                       * zones
                                       */
Pout<< "zoneName:" << zoneName
    //<< " varLocArray:" << varLocArray
    << " solTime:" << SolTime
    << endl;



    INTEGER4 *PassiveVarArray = nullptr;
    INTEGER4 *VarShareArray   = nullptr;
    INTEGER4  ShrConn         = 0;

    INTEGER4  NumBConns       = 0;   /* No Boundary Connections */
    INTEGER4  NumBItems       = 0;   /* No Boundary Items */

    INTEGER4  ZoneType = ZoneType_FEPolyhedron;

    if
    (
       !TECZNE112
        (
            const_cast<char*>(zoneName.c_str()),
            &ZoneType,
            &NumNodes,
            &NumElems,
            &NumFaces,
            &ICellMax,
            &JCellMax,
            &KCellMax,
            &SolTime,
            &strandID,
            &ParentZone,
            &IsBlock,
            &NFConns,
            &FNMode,
            &nFaceNodes,
            &NumBConns,
            &NumBItems,
            PassiveVarArray,
            const_cast<INTEGER4*>(varLocArray.begin()),
            VarShareArray,
            &ShrConn
        )
    )
    {
//        FatalErrorInFunction
//            << "Error in TECZNE112." << exit(FatalError);
    }
}


void Foam::tecplotWriter::writePolygonalZone
(
    const word& zoneName,
    INTEGER4 strandID,
    const indirectPrimitivePatch& pp,
    const List<INTEGER4>& varLocArray
) const
{
    /* Call TECZNE112 */
    INTEGER4  NumNodes   = pp.nPoints();         /* number of unique nodes */
    INTEGER4  NumElems   = pp.size();         /* number of elements */
    INTEGER4  NumFaces   = pp.nEdges();         /* number of unique faces */

    INTEGER4  ICellMax   = 0;         /* Not Used, set to zero */
    INTEGER4  JCellMax   = 0;         /* Not Used, set to zero */
    INTEGER4  KCellMax   = 0;         /* Not Used, set to zero */

    double    SolTime    = runTime_.value();     /* solution time   */
    INTEGER4  ParentZone = 0;         /* no parent zone  */

    INTEGER4  IsBlock    = 1;         /* block format  */

    INTEGER4  NFConns    = 0;         /* not used for FEPolyhedron
                                       * zones
                                       */
    INTEGER4  FNMode     = 0;         /* not used for FEPolyhedron
                                       * zones
                                       */
    INTEGER4  NumFaceNodes    = 2*pp.nEdges();

Pout<< "zoneName:" << zoneName
    << " strandID:" << strandID
    //<< " varLocArray:" << varLocArray
    << " solTime:" << SolTime
    << endl;


    INTEGER4 *PassiveVarArray = nullptr;
    INTEGER4 *VarShareArray   = nullptr;
    INTEGER4  ShrConn         = 0;

    INTEGER4  NumBConns       = 0;   /* No Boundary Connections */
    INTEGER4  NumBItems       = 0;   /* No Boundary Items */

    INTEGER4  ZoneType = ZoneType_FEPolygon;

    if
    (
       !TECZNE112
        (
            const_cast<char*>(zoneName.c_str()),
            &ZoneType,
            &NumNodes,
            &NumElems,
            &NumFaces,
            &ICellMax,
            &JCellMax,
            &KCellMax,
            &SolTime,
            &strandID,
            &ParentZone,
            &IsBlock,
            &NFConns,
            &FNMode,
            &NumFaceNodes,
            &NumBConns,
            &NumBItems,
            PassiveVarArray,
            const_cast<INTEGER4*>(varLocArray.begin()),
            VarShareArray,
            &ShrConn
        )
    )
    {
//        FatalErrorInFunction
//            << "Error in TECZNE112." << exit(FatalError);
    }
}


void Foam::tecplotWriter::writeOrderedZone
(
    const word& zoneName,
    INTEGER4 strandID,
    const label n,
    const List<INTEGER4>& varLocArray
) const
{
    /* Call TECZNE112 */
    INTEGER4  IMax   = n;         /* number of unique nodes */
    INTEGER4  JMax   = 1;         /* number of elements */
    INTEGER4  KMax   = 1;         /* number of unique faces */

    INTEGER4  ICellMax   = 0;         /* Not Used, set to zero */
    INTEGER4  JCellMax   = 0;         /* Not Used, set to zero */
    INTEGER4  KCellMax   = 0;         /* Not Used, set to zero */

    double    SolTime    = runTime_.value();     /* solution time   */
    INTEGER4  ParentZone = 0;         /* no parent zone  */

    INTEGER4  IsBlock    = 1;         /* block format  */

    INTEGER4  NFConns    = 0;         /* not used for FEPolyhedron
                                       * zones
                                       */
    INTEGER4  FNMode     = 0;         /* not used for FEPolyhedron
                                       * zones
                                       */
    INTEGER4  NumFaceNodes    = 1;
    INTEGER4  NumBConns       = 1;   /* No Boundary Connections */
    INTEGER4  NumBItems       = 1;   /* No Boundary Items */

Pout<< "zoneName:" << zoneName
    << " strandID:" << strandID
    //<< " varLocArray:" << varLocArray
    << " solTime:" << SolTime
    << endl;


    INTEGER4 *PassiveVarArray = nullptr;
    INTEGER4 *VarShareArray   = nullptr;
    INTEGER4  ShrConn         = 0;


    INTEGER4  ZoneType = ZoneType_Ordered;

    if
    (
       !TECZNE112
        (
            const_cast<char*>(zoneName.c_str()),
            &ZoneType,
            &IMax,
            &JMax,
            &KMax,
            &ICellMax,
            &JCellMax,
            &KCellMax,
            &SolTime,
            &strandID,
            &ParentZone,
            &IsBlock,
            &NFConns,
            &FNMode,
            &NumFaceNodes,
            &NumBConns,
            &NumBItems,
            PassiveVarArray,
            const_cast<INTEGER4*>(varLocArray.begin()),
            VarShareArray,
            &ShrConn
        )
    )
    {
//        FatalErrorInFunction
//            << "Error in TECZNE112." << exit(FatalError);
    }
}


void Foam::tecplotWriter::writeConnectivity(const fvMesh& mesh) const
{
    List<INTEGER4> FaceNodeCounts(mesh.nFaces());

    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];
        FaceNodeCounts[facei] = INTEGER4(f.size());
    }


    INTEGER4 nFaceNodes = 0;
    forAll(mesh.faces(), facei)
    {
        nFaceNodes += mesh.faces()[facei].size();
    }


    List<INTEGER4> FaceNodes(nFaceNodes);
    label nodeI = 0;
    forAll(mesh.faces(), facei)
    {
        const face& f = mesh.faces()[facei];
        forAll(f, fp)
        {
            FaceNodes[nodeI++] = INTEGER4(f[fp]+1);
        }
    }


    List<INTEGER4> FaceLeftElems(mesh.nFaces());
    forAll(mesh.faceOwner(), facei)
    {
        FaceLeftElems[facei] = mesh.faceOwner()[facei]+1;
    }

    List<INTEGER4> FaceRightElems(mesh.nFaces());
    forAll(mesh.faceNeighbour(), facei)
    {
        FaceRightElems[facei] = mesh.faceNeighbour()[facei]+1;
    }
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        FaceRightElems[facei] = 0;
    }

    if
    (
       !TECPOLY112
        (
            FaceNodeCounts.begin(), /* The face node counts array */
            FaceNodes.begin(),      /* The face nodes array */
            FaceLeftElems.begin(),  /* The left elements array  */
            FaceRightElems.begin(), /* The right elements array  */
            nullptr,       /* No boundary connection counts */
            nullptr,       /* No boundary connection elements */
            nullptr        /* No boundary connection zones */
        )
    )
    {
//        FatalErrorInFunction
//            << "Error in TECPOLY112." << exit(FatalError);
    }
}

void Foam::tecplotWriter::writeConnectivity
(
    const indirectPrimitivePatch& pp
) const
{
    INTEGER4  NumFaces   = pp.nEdges();         /* number of unique faces */
    INTEGER4  NumFaceNodes    = 2*pp.nEdges();

    // All faces (=edges) have 2 nodes
    List<INTEGER4> FaceNodeCounts(NumFaces);
    FaceNodeCounts = 2;

    List<INTEGER4> FaceNodes(NumFaceNodes);
    label nodeI = 0;
    forAll(pp.edges(), edgeI)
    {
        edge e = pp.edges()[edgeI];
        if (e[0] > e[1])
        {
            e.flip();
        }

        FaceNodes[nodeI++] = INTEGER4(e[0]+1);
        FaceNodes[nodeI++] = INTEGER4(e[1]+1);
    }

    /* Define the right and left elements of each face.
     *
     * The last step for writing out the polyhedral data is to
     * define the right and left neighboring elements for each
     * face.  The neighboring elements can be determined using the
     * right-hand rule.  For each face, place your right-hand along
     * the face which your fingers pointing the direction of
     * incrementing node numbers (i.e. from node 1 to node 2).
     * Your right thumb will point towards the right element; the
     * element on the other side of your hand is the left element.
     *
     * The number zero is used to indicate that there isn't an
     * element on that side of the face.
     *
     * Because of the way we numbered the nodes and faces, the
     * right element for every face is the element itself
     * (element 1) and the left element is "no-neighboring element"
     * (element 0).
     */

    List<INTEGER4> FaceLeftElems(NumFaces);
    List<INTEGER4> FaceRightElems(NumFaces);

    const labelListList& edgeFaces = pp.edgeFaces();
    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 1)
        {
            FaceLeftElems[edgeI] = 0;
            FaceRightElems[edgeI] = eFaces[0]+1;
        }
        else if (eFaces.size() == 2)
        {
            edge e = pp.edges()[edgeI];
            if (e[0] > e[1])
            {
                e.flip();
            }

            const face& f0 = pp.localFaces()[eFaces[0]];

            // The face that uses the vertices of e in increasing order
            // is the left face.

            label fp = findIndex(f0, e[0]);
            bool f0IsLeft = (f0.nextLabel(fp) == e[1]);

            if (f0IsLeft)
            {
                FaceLeftElems[edgeI] = eFaces[0]+1;
                FaceRightElems[edgeI] = eFaces[1]+1;
            }
            else
            {
                FaceLeftElems[edgeI] = eFaces[1]+1;
                FaceRightElems[edgeI] = eFaces[0]+1;
            }
        }
        else
        {
            // non-manifold. Treat as if open.
            FaceLeftElems[edgeI] = 0;
            FaceRightElems[edgeI] = eFaces[0]+1;
        }
    }

    /* Write the face map (created above) using TECPOLY112. */
    if
    (
       !TECPOLY112
        (
            FaceNodeCounts.begin(), /* The face node counts array */
            FaceNodes.begin(),      /* The face nodes array */
            FaceLeftElems.begin(),  /* The left elements array  */
            FaceRightElems.begin(), /* The right elements array  */
            nullptr,       /* No boundary connection counts */
            nullptr,       /* No boundary connection elements */
            nullptr        /* No boundary connection zones */
        )
    )
    {
//        FatalErrorInFunction
//            << "Error in TECPOLY112." << exit(FatalError);
    }
}


void Foam::tecplotWriter::writeEnd() const
{
Pout<< "writeEnd" << endl;

    if (!TECEND112())
    {
//        FatalErrorInFunction
//            << "Error in TECEND112." << exit(FatalError);
    }

}


// ************************************************************************* //
