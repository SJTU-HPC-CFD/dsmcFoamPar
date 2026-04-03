/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2506                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
Build  : 95a5dfacc7-20250714 OPENFOAM=2506 version=v2506
Arch   : "LSB;label=32;scalar=64"
Exec   : /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/platforms/linux64IcxDPInt32Opt/bin/dsmcFoam+
Date   : Apr 02 2026
Time   : 23:54:46
Host   : SuperXCX-ROG
PID    : 2398822
I/O    : uncollated
Case   : /home/superxcx/code/DSMC/dsmcFoam++/hyStrath_xcx/case/cylinderN2noreact_test
nProcs : 1
trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 5, maxFileModificationPolls 20)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

--> FOAM Warning : 
    From static Foam::IOstreamOption::compressionType Foam::IOstreamOption::compressionEnum(const word &, const compressionType)
    in file db/options/IOstreamOption.C at line 189
    Unknown compression specifier 'uncompressed' using compression off
Create mesh for time = 0


Constructing dsmcCloud
Selecting the coordinate system model:	dsmcCartesian

Selecting the time-step model:	dsmcConstantTimeStepModel

Constant time-step model:
- time-step [sec]	6.632426865e-08

Selecting the porous measurement model:	dsmcNoPorousMediumMeasurements


Creating dsmcControllers

Creating fields: 

Selecting field: dsmcVolFields

TimeData Statistics: 
Resetting fields at output until time: 4.841671611e-06
 measurement option: write
 nSamples: 1, time interval: 6.632426865e-08
 nAverages: 83, time interval: 5.504914298e-06
 total no. of sampling steps: 0
 total no. of averaging Steps: 0


Creating the boundary models: 

Selecting dsmcPatchBoundaryModel dsmcDeletionPatch
Selecting dsmcPatchBoundaryModel dsmcDeletionPatch
Selecting dsmcPatchBoundaryModel dsmcDiffuseWallPatch
Selecting dsmcGeneralBoundaryModel dsmcFreeStreamInflowPatch
Selecting BinaryCollisionModel LarsenBorgnakkeVariableHardSphere
OpenMP enabled for dsmcCloud with 1 thread-local RNG streams using collision strategy 'partition'

Constructing constant properties for
    N2

Creating dsmcReactions

There are no chemical reactions defined.
Selecting collisionPartnerSelectionModel noTimeCounter

Initialising the measurement fields

Initialising dsmcVolFields field

Starting time loop

ExecutionTime = 10.1 s  ClockTime = 9 s

ExecutionTime = 14.18 s  ClockTime = 11 s

ExecutionTime = 18.43 s  ClockTime = 10 s

ExecutionTime = 22.48 s  ClockTime = 11 s

ExecutionTime = 26.7 s  ClockTime = 13 s

ExecutionTime = 30.9 s  ClockTime = 15 s

ExecutionTime = 34.99 s  ClockTime = 17 s

ExecutionTime = 39.21 s  ClockTime = 18 s

ExecutionTime = 43.62 s  ClockTime = 20 s

Time = 6.632426865e-07

    Collisions                      = 28198

Cloud name: dsmc
    Number of dsmc particles        = 1329896
    Number of molecules             = 7.979376e+15
    Mass in system                  = 3.71040984e-10
    Average linear kinetic energy   = 1.836409341e-19
    Average rotational energy       = 1.338316505e-21
    Average vibrational energy      = 2.327089075e-20
    Average electronic energy       = 0
    Average total energy            = 2.082501413e-19
ExecutionTime = 47.95 s  ClockTime = 22 s

ExecutionTime = 52.07 s  ClockTime = 24 s

ExecutionTime = 56.4 s  ClockTime = 26 s

ExecutionTime = 60.77 s  ClockTime = 28 s

ExecutionTime = 65.02 s  ClockTime = 30 s

ExecutionTime = 69.43 s  ClockTime = 31 s

ExecutionTime = 73.65 s  ClockTime = 33 s

ExecutionTime = 78.11 s  ClockTime = 35 s

ExecutionTime = 82.41 s  ClockTime = 37 s

ExecutionTime = 86.8 s  ClockTime = 39 s

Time = 1.326485373e-06

    Collisions                      = 37431

Cloud name: dsmc
    Number of dsmc particles        = 1353515
    Number of molecules             = 8.12109e+15
    Mass in system                  = 3.77630685e-10
    Average linear kinetic energy   = 1.809908722e-19
    Average rotational energy       = 1.572197159e-21
    Average vibrational energy      = 2.327085575e-20
    Average electronic energy       = 0
    Average total energy            = 2.058339251e-19
ExecutionTime = 91.53 s  ClockTime = 38 s

ExecutionTime = 95.83 s  ClockTime = 41 s

ExecutionTime = 100.31 s  ClockTime = 43 s

ExecutionTime = 104.7 s  ClockTime = 45 s

ExecutionTime = 109.19 s  ClockTime = 47 s

ExecutionTime = 113.68 s  ClockTime = 49 s

ExecutionTime = 118.24 s  ClockTime = 51 s

ExecutionTime = 122.78 s  ClockTime = 53 s

ExecutionTime = 127.24 s  ClockTime = 55 s

ExecutionTime = 131.92 s  ClockTime = 57 s

Time = 1.98972806e-06

    Collisions                      = 47932

Cloud name: dsmc
    Number of dsmc particles        = 1377011
    Number of molecules             = 8.262066e+15
    Mass in system                  = 3.84186069e-10
    Average linear kinetic energy   = 1.786399273e-19
    Average rotational energy       = 1.889420419e-21
    Average vibrational energy      = 2.327085575e-20
    Average electronic energy       = 0
    Average total energy            = 2.038002035e-19
ExecutionTime = 136.64 s  ClockTime = 59 s

ExecutionTime = 140.94 s  ClockTime = 61 s

ExecutionTime = 145.74 s  ClockTime = 64 s

ExecutionTime = 150.61 s  ClockTime = 66 s

ExecutionTime = 155.26 s  ClockTime = 69 s

ExecutionTime = 160.01 s  ClockTime = 68 s

ExecutionTime = 164.59 s  ClockTime = 70 s

ExecutionTime = 169.41 s  ClockTime = 73 s

ExecutionTime = 174.16 s  ClockTime = 75 s

ExecutionTime = 178.82 s  ClockTime = 77 s

Time = 2.652970746e-06

    Collisions                      = 58096

Cloud name: dsmc
    Number of dsmc particles        = 1400432
    Number of molecules             = 8.402592e+15
    Mass in system                  = 3.90720528e-10
    Average linear kinetic energy   = 1.765076833e-19
    Average rotational energy       = 2.285665764e-21
    Average vibrational energy      = 2.327092222e-20
    Average electronic energy       = 0
    Average total energy            = 2.020642713e-19
ExecutionTime = 183.94 s  ClockTime = 80 s

ExecutionTime = 188.42 s  ClockTime = 82 s

ExecutionTime = 193.14 s  ClockTime = 85 s

ExecutionTime = 197.9 s  ClockTime = 87 s

ExecutionTime = 202.5 s  ClockTime = 89 s

ExecutionTime = 207.3 s  ClockTime = 92 s

ExecutionTime = 212.11 s  ClockTime = 94 s

ExecutionTime = 216.93 s  ClockTime = 97 s

ExecutionTime = 221.65 s  ClockTime = 99 s

ExecutionTime = 226.69 s  ClockTime = 99 s

Time = 3.316213433e-06

    Collisions                      = 68592

Cloud name: dsmc
    Number of dsmc particles        = 1423886
    Number of molecules             = 8.543316e+15
    Mass in system                  = 3.97264194e-10
    Average linear kinetic energy   = 1.745672564e-19
    Average rotational energy       = 2.724507917e-21
    Average vibrational energy      = 2.32709865e-20
    Average electronic energy       = 0
    Average total energy            = 2.005627508e-19
ExecutionTime = 231.58 s  ClockTime = 101 s

ExecutionTime = 236.23 s  ClockTime = 104 s

ExecutionTime = 241.24 s  ClockTime = 106 s

ExecutionTime = 246.04 s  ClockTime = 109 s

ExecutionTime = 250.91 s  ClockTime = 111 s

ExecutionTime = 255.97 s  ClockTime = 114 s

ExecutionTime = 260.91 s  ClockTime = 116 s

ExecutionTime = 265.87 s  ClockTime = 119 s

ExecutionTime = 270.73 s  ClockTime = 122 s

ExecutionTime = 275.84 s  ClockTime = 124 s

Time = 3.979456119e-06

    Collisions                      = 77899

Cloud name: dsmc
    Number of dsmc particles        = 1447268
    Number of molecules             = 8.683608e+15
    Mass in system                  = 4.03787772e-10
    Average linear kinetic energy   = 1.72772351e-19
    Average rotational energy       = 3.188737097e-21
    Average vibrational energy      = 2.32710487e-20
    Average electronic energy       = 0
    Average total energy            = 1.992321368e-19
ExecutionTime = 280.86 s  ClockTime = 127 s

ExecutionTime = 285.76 s  ClockTime = 129 s

ExecutionTime = 290.76 s  ClockTime = 129 s

ExecutionTime = 295.72 s  ClockTime = 132 s

ExecutionTime = 300.77 s  ClockTime = 134 s

ExecutionTime = 305.82 s  ClockTime = 137 s

ExecutionTime = 310.86 s  ClockTime = 140 s

ExecutionTime = 315.92 s  ClockTime = 143 s

ExecutionTime = 320.97 s  ClockTime = 145 s

ExecutionTime = 326.01 s  ClockTime = 148 s

Time = 4.642698805e-06

    Collisions                      = 87829

Cloud name: dsmc
    Number of dsmc particles        = 1470593
    Number of molecules             = 8.823558e+15
    Mass in system                  = 4.10295447e-10
    Average linear kinetic energy   = 1.711041089e-19
    Average rotational energy       = 3.669968608e-21
    Average vibrational energy      = 2.327142542e-20
    Average electronic energy       = 0
    Average total energy            = 1.980455029e-19
ExecutionTime = 331.17 s  ClockTime = 151 s

ExecutionTime = 336.24 s  ClockTime = 154 s


Start sample at time 4.775347343e-06 and step 72

ExecutionTime = 341.33 s  ClockTime = 156 s

ExecutionTime = 346.42 s  ClockTime = 159 s

ExecutionTime = 351.51 s  ClockTime = 159 s

ExecutionTime = 356.62 s  ClockTime = 162 s

ExecutionTime = 361.84 s  ClockTime = 165 s

ExecutionTime = 366.98 s  ClockTime = 167 s

ExecutionTime = 372.15 s  ClockTime = 170 s

ExecutionTime = 377.3 s  ClockTime = 173 s

Time = 5.305941492e-06

    Collisions                      = 96787

Cloud name: dsmc
    Number of dsmc particles        = 1494103
    Number of molecules             = 8.964618e+15
    Mass in system                  = 4.16854737e-10
    Average linear kinetic energy   = 1.695323787e-19
    Average rotational energy       = 4.178133946e-21
    Average vibrational energy      = 2.327129185e-20
    Average electronic energy       = 0
    Average total energy            = 1.969818045e-19
ExecutionTime = 382.54 s  ClockTime = 176 s

ExecutionTime = 387.74 s  ClockTime = 179 s

ExecutionTime = 392.92 s  ClockTime = 182 s

BuildCellOccupancy profiling summary:
    buildCellOccupancy calls      = 84
    extract parcels [s]           = 4.728131306
    count/reduce [s]              = 0
    allocate/fill [s]             = 0
    total profiled [s]            = 4.728131306

Collision profiling summary:
    collision calls               = 83
    precompute candidates [s]     = 0.460595276
    rebuild partition [s]         = 6.8655e-05
    selection/collide [s]         = 57.38651911
    total profiled [s]            = 57.84718304


Sample average steps = 11

Total wall force: 
 wallForce_x 0.02375831899
 wallForce_y -0.006898064908
 wallForce_z -5.303055196e-05
Evolve profiling summary:
    evolve calls                  = 83
    move only [s]                 = 117.9729657
    buildCellOccupancy [s]        = 4.701674454
    coordSystem [s]               = 0.000206507
    collision phase [s]           = 57.84817329
    reaction/output [s]           = 1.5268e-05
    post fields/output [s]        = 14.05765146
    total profiled [s]            = 194.5806867

ExecutionTime = 412.94 s  ClockTime = 198 s

End

