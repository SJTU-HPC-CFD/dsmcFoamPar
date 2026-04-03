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
Time   : 20:47:57
Host   : SuperXCX-ROG
PID    : 2363840
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
OpenMP enabled for dsmcCloud with 8 thread-local RNG streams

Constructing constant properties for
    N2

Creating dsmcReactions

There are no chemical reactions defined.
Selecting collisionPartnerSelectionModel noTimeCounter

Initialising the measurement fields

Initialising dsmcVolFields field

Starting time loop

ExecutionTime = 10.16 s  ClockTime = 10 s

ExecutionTime = 13.17 s  ClockTime = 11 s

ExecutionTime = 16.21 s  ClockTime = 13 s

ExecutionTime = 19.33 s  ClockTime = 14 s

ExecutionTime = 22.47 s  ClockTime = 16 s

ExecutionTime = 25.52 s  ClockTime = 17 s

ExecutionTime = 28.84 s  ClockTime = 19 s

ExecutionTime = 32.08 s  ClockTime = 20 s

ExecutionTime = 35.33 s  ClockTime = 22 s

Time = 6.632426865e-07

    Collisions                      = 28195

Cloud name: dsmc
    Number of dsmc particles        = 1329948
    Number of molecules             = 7.979688e+15
    Mass in system                  = 3.71055492e-10
    Average linear kinetic energy   = 1.836289486e-19
    Average rotational energy       = 1.338106817e-21
    Average vibrational energy      = 2.327089075e-20
    Average electronic energy       = 0
    Average total energy            = 2.082379462e-19
ExecutionTime = 38.68 s  ClockTime = 23 s

ExecutionTime = 41.64 s  ClockTime = 22 s

ExecutionTime = 44.8 s  ClockTime = 24 s

ExecutionTime = 47.96 s  ClockTime = 25 s

ExecutionTime = 51.34 s  ClockTime = 27 s

ExecutionTime = 54.64 s  ClockTime = 28 s

ExecutionTime = 57.79 s  ClockTime = 30 s

ExecutionTime = 61.37 s  ClockTime = 31 s

ExecutionTime = 64.78 s  ClockTime = 33 s

ExecutionTime = 68.08 s  ClockTime = 34 s

Time = 1.326485373e-06

    Collisions                      = 37253

Cloud name: dsmc
    Number of dsmc particles        = 1353437
    Number of molecules             = 8.120622e+15
    Mass in system                  = 3.77608923e-10
    Average linear kinetic energy   = 1.809910978e-19
    Average rotational energy       = 1.573702507e-21
    Average vibrational energy      = 2.327085575e-20
    Average electronic energy       = 0
    Average total energy            = 2.05835656e-19
ExecutionTime = 72.09 s  ClockTime = 36 s

ExecutionTime = 75.08 s  ClockTime = 38 s

ExecutionTime = 78.63 s  ClockTime = 39 s

ExecutionTime = 82.21 s  ClockTime = 41 s

ExecutionTime = 85.81 s  ClockTime = 43 s

ExecutionTime = 89.34 s  ClockTime = 45 s

ExecutionTime = 93.26 s  ClockTime = 47 s

ExecutionTime = 97.01 s  ClockTime = 48 s

ExecutionTime = 100.91 s  ClockTime = 50 s

ExecutionTime = 104.61 s  ClockTime = 52 s

Time = 1.98972806e-06

    Collisions                      = 48155

Cloud name: dsmc
    Number of dsmc particles        = 1376995
    Number of molecules             = 8.26197e+15
    Mass in system                  = 3.84181605e-10
    Average linear kinetic energy   = 1.786444491e-19
    Average rotational energy       = 1.896717566e-21
    Average vibrational energy      = 2.327092335e-20
    Average electronic energy       = 0
    Average total energy            = 2.038120901e-19
ExecutionTime = 108.49 s  ClockTime = 54 s

ExecutionTime = 111.93 s  ClockTime = 53 s

ExecutionTime = 115.63 s  ClockTime = 55 s

ExecutionTime = 119.38 s  ClockTime = 57 s

ExecutionTime = 123.12 s  ClockTime = 58 s

ExecutionTime = 126.85 s  ClockTime = 60 s

ExecutionTime = 130.7 s  ClockTime = 62 s

ExecutionTime = 134.57 s  ClockTime = 64 s

ExecutionTime = 138.36 s  ClockTime = 66 s

ExecutionTime = 142.28 s  ClockTime = 68 s

Time = 2.652970746e-06

    Collisions                      = 57566

Cloud name: dsmc
    Number of dsmc particles        = 1400530
    Number of molecules             = 8.40318e+15
    Mass in system                  = 3.9074787e-10
    Average linear kinetic energy   = 1.765086682e-19
    Average rotational energy       = 2.285206037e-21
    Average vibrational energy      = 2.327095545e-20
    Average electronic energy       = 0
    Average total energy            = 2.020648297e-19
ExecutionTime = 146.59 s  ClockTime = 70 s

ExecutionTime = 150.13 s  ClockTime = 72 s

ExecutionTime = 154.17 s  ClockTime = 74 s

ExecutionTime = 158.01 s  ClockTime = 76 s

ExecutionTime = 161.98 s  ClockTime = 77 s

ExecutionTime = 166.05 s  ClockTime = 79 s

ExecutionTime = 170.33 s  ClockTime = 82 s

ExecutionTime = 174.26 s  ClockTime = 83 s

ExecutionTime = 178.29 s  ClockTime = 83 s

ExecutionTime = 182.38 s  ClockTime = 85 s

Time = 3.316213433e-06

    Collisions                      = 68294

Cloud name: dsmc
    Number of dsmc particles        = 1424005
    Number of molecules             = 8.54403e+15
    Mass in system                  = 3.97297395e-10
    Average linear kinetic energy   = 1.745618676e-19
    Average rotational energy       = 2.724062881e-21
    Average vibrational energy      = 2.32709538e-20
    Average electronic energy       = 0
    Average total energy            = 2.005568843e-19
ExecutionTime = 186.82 s  ClockTime = 87 s

ExecutionTime = 190.67 s  ClockTime = 89 s

ExecutionTime = 194.7 s  ClockTime = 91 s

ExecutionTime = 198.85 s  ClockTime = 93 s

ExecutionTime = 203.2 s  ClockTime = 95 s

ExecutionTime = 207.56 s  ClockTime = 97 s

ExecutionTime = 211.59 s  ClockTime = 99 s

ExecutionTime = 216.01 s  ClockTime = 101 s

ExecutionTime = 220.54 s  ClockTime = 103 s

ExecutionTime = 224.99 s  ClockTime = 105 s

Time = 3.979456119e-06

    Collisions                      = 78099

Cloud name: dsmc
    Number of dsmc particles        = 1447292
    Number of molecules             = 8.683752e+15
    Mass in system                  = 4.03794468e-10
    Average linear kinetic energy   = 1.727668931e-19
    Average rotational energy       = 3.180994741e-21
    Average vibrational energy      = 2.327095222e-20
    Average electronic energy       = 0
    Average total energy            = 1.9921884e-19
ExecutionTime = 229.55 s  ClockTime = 108 s

ExecutionTime = 233.75 s  ClockTime = 110 s

ExecutionTime = 237.92 s  ClockTime = 112 s

ExecutionTime = 242.21 s  ClockTime = 114 s

ExecutionTime = 246.45 s  ClockTime = 113 s

ExecutionTime = 250.98 s  ClockTime = 115 s

ExecutionTime = 255.19 s  ClockTime = 118 s

ExecutionTime = 260.05 s  ClockTime = 120 s

ExecutionTime = 264.34 s  ClockTime = 122 s

ExecutionTime = 268.7 s  ClockTime = 124 s

Time = 4.642698805e-06

    Collisions                      = 87633

Cloud name: dsmc
    Number of dsmc particles        = 1470734
    Number of molecules             = 8.824404e+15
    Mass in system                  = 4.10334786e-10
    Average linear kinetic energy   = 1.710950614e-19
    Average rotational energy       = 3.672093543e-21
    Average vibrational energy      = 2.327107727e-20
    Average electronic energy       = 0
    Average total energy            = 1.980382322e-19
ExecutionTime = 273.45 s  ClockTime = 127 s

ExecutionTime = 277.95 s  ClockTime = 129 s


Start sample at time 4.775347343e-06 and step 72

ExecutionTime = 282.52 s  ClockTime = 131 s

ExecutionTime = 286.86 s  ClockTime = 133 s

ExecutionTime = 291.32 s  ClockTime = 136 s

ExecutionTime = 295.82 s  ClockTime = 138 s

ExecutionTime = 300.51 s  ClockTime = 140 s

ExecutionTime = 305 s  ClockTime = 143 s

ExecutionTime = 309.61 s  ClockTime = 142 s

ExecutionTime = 314.11 s  ClockTime = 144 s

Time = 5.305941492e-06

    Collisions                      = 95650

Cloud name: dsmc
    Number of dsmc particles        = 1494110
    Number of molecules             = 8.96466e+15
    Mass in system                  = 4.1685669e-10
    Average linear kinetic energy   = 1.69510249e-19
    Average rotational energy       = 4.176498365e-21
    Average vibrational energy      = 2.32713853e-20
    Average electronic energy       = 0
    Average total energy            = 1.969581327e-19
ExecutionTime = 318.61 s  ClockTime = 147 s

ExecutionTime = 323.21 s  ClockTime = 149 s

ExecutionTime = 327.79 s  ClockTime = 151 s

Collision profiling summary:
    collision calls               = 83
    precompute candidates [s]     = 0.159875935
    rebuild partition [s]         = 0.012676563
    selection/collide [s]         = 10.96253542
    total profiled [s]            = 11.13508791


Sample average steps = 11

Total wall force: 
 wallForce_x 0.0235693621
 wallForce_y -0.006957683365
 wallForce_z -1.92653176e-05
ExecutionTime = 347.71 s  ClockTime = 170 s

End

