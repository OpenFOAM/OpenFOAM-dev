const pointField& points = mesh.points();
label nPoints = points.size();
const cellShapeList& cells = mesh.cellShapes();

gmvFile << "gmvinput " << format << nl;
gmvFile << "nodes   " << nPoints << nl;
for (label indx=0;indx<nPoints;indx++)
{
    gmvFile << points[indx].x() << " ";
}
gmvFile << nl;
for (label indx=0;indx<nPoints;indx++)
{
    gmvFile << points[indx].y() << " ";
}
gmvFile << nl;
for (label indx=0;indx<nPoints;indx++)
{
    gmvFile << points[indx].z() << " ";
}
gmvFile << nl;
gmvFile << "cells   " << cells.size() << nl;
forAll(cells, indx)
{
    label nNodes = cells[indx].size();
    if (nNodes == 8)
    {
        gmvFile << "hex     " << 8 << " ";
        for (label ip=0; ip<nNodes; ip++)
        {
            gmvFile << cells[indx][ip] + 1 << " ";
        }
        gmvFile << nl;
    }
    else if (nNodes == 4)
    {
        gmvFile << "tet     " << 4 << " ";
        for (label ip=0; ip<nNodes; ip++)
        {
            gmvFile << cells[indx][ip] + 1 << " ";
        }
        gmvFile << nl;
    }
    else if (nNodes == 6)
    {
        gmvFile << "prism   " << 6 << " ";
        for (label ip=0; ip<nNodes; ip++)
        {
            gmvFile << cells[indx][ip] + 1 << " ";
        }
        gmvFile << nl;
    }
}
