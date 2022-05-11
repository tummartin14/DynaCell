;--------INPUT--------
dc_modus  md
pdbfile  protein.pdb
ligandpdbfile  md/lig.pdb               ;if ligand-pdb is a separate file
topologyfile  protein.top
ligandtopologyfile  ligand.top          ;if ligand-topology is a separate file
backboneFile  md/backboneNames.dat      ;needed if softcore/restraints need bb/sc definiton

;--------OUTPUT--------
outputfile  md/md.log.out
trajectoryfile  md/traj.out.dcd
trajwrite  400                          ;frequency of trajectory output
outwrite  400                           ;frequency of log output
nrsteps  60000                          ;maximum nr of steps if SC-alpha does not reach zero
randomseed  1234567

;--------RESTRAINTS--------
restraintAtomsOuter  backbone           ;restraints to be set in the outer region
restraintAtomsInner  none               ;restraints to be set in the inner region
restraintForceconstOuter  1000.0        ;Fc for outer region, in kJ/(mol*nm^2)
restraintRegionLigandRadius  1.2        ;define inner region as ligand and the surrounding area, in nm

;--------ENERGY--------
solvent  AGBNP2
AGBNPparaFile  md/AGBNP2.dat
solventDielectricConst  80.0
solventScaling  excludeLig              ;enable solvent-scaling feature by
solventScalingAlphaUp  0.0              ;coupling the solvent model evaluation
solventScalingAlphaLow  0.0             ;to the SC-alpha value
useChargeGroups  no

;--------MD--------
langevin  yes                           ;use Langevin dynamics
frictionCoefficient  0.5                ;in 1/ps
temperature  300                        ;in Kelvin
timestep  0.001                         ;in ps
constraints  shake
constraintsType  Hbonds
VcomRemoval  100                        ;center of mass removal frequency
nrStepsAfterSCisDone  40000             ;nr of steps to simulate after SC-alpha is zero
resetOutWritesAfterSCisDone  yes
restartfileout  md/restart.out

;--------SOFTCORE--------
SCprogression  optimize                 ;determine SC-alpha via minimization 
SCprogressionPots  LJ_Q                 ;set SC on non-bonded potentials
SCnbRegion  ligandIntAct                ;apply SC-potentials between ligand and protein
SCinitialAlphaLJ  -1.0
SCinitialAlphaQ  -1.0
SCautomaticFrequency  yes
SCOptCoupleAlphas  yes
