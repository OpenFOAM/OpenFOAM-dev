//======================================================================
// Global variables
const int maxNames = 1000;

// define the name to be displayed in the window.
static char readerName[]           = "OpenFOAM";
static char meshName[]             = "cells";
static char readerVersion[]        = "2.03";

// everything is one part in OpenFOAM, except the spray
static int Num_unstructured_parts  = 1;
static int Num_structured_parts    = 0;
static int Numparts_available      = 1;
static int nPatches                = 0;

static int Num_timesets            = 1;
static int Geom_timeset_number     = 1;
static int Num_time_steps          = 1;
static int Num_global_nodes        = 0;
static int Num_variables           = 0;
static int Num_dataset_files       = 0;
static int Current_time_step       = 0;

static label nSprayVariables       = 0;
static label nMaxParcels           = 0;

static bool isScalar[maxNames];
static bool isVector[maxNames];
static bool isTensor[maxNames];
static bool isSpray[maxNames];

static word scalarName             = "volScalarField";
static word vectorName             = "volVectorField";
static word tensorName             = "volTensorField";
static word sprayScalarFieldName   = "scalarField";
static word sprayVectorFieldName   = "vectorField";
static word sprayTensorFieldName   = "tensorField";
static word parcelPrepend          = "parcel_";
static word pointPrepend           = "point_";

static fileName rootDir;
static fileName caseDir;

static instantList timeDirs;

static List<word> fieldNames;
static List<word> lagrangianScalarNames;
static List<word> lagrangianVectorNames;
static label var2field[maxNames];

static Time *runTimePtr = 0;
static fvMesh *meshPtr = 0;
static Cloud<passiveParticle> *sprayPtr = 0;
