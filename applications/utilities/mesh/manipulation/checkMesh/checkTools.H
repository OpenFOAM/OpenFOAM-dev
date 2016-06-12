#include "scalar.H"
#include "indirectPrimitivePatch.H"

namespace Foam
{
    class polyMesh;
    class surfaceWriter;
    class faceSet;
    class cellSet;
    class fileName;
    class polyMesh;

    void printMeshStats(const polyMesh& mesh, const bool allTopology);

    //- Generate merged surface on master and write. Needs input patch
    //  to be of mesh faces.
    void mergeAndWrite
    (
        const polyMesh& mesh,
        const surfaceWriter& writer,
        const word& name,
        const indirectPrimitivePatch setPatch,
        const fileName& outputDir
    );

    //- Write vtk representation of (assembled) faceSet to vtk file in
    //  postProcessing/ directory
    void mergeAndWrite(const surfaceWriter&, const faceSet&);

    //- Write vtk representation of (assembled) cellSet to vtk file in
    //  postProcessing/ directory
    void mergeAndWrite(const surfaceWriter&, const cellSet&);
}
