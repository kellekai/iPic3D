Module Constants
    Implicit none
    Real(8),parameter :: PI=3.141592653589793238D0
    Real(8),parameter :: a0=5.2918d-11
    Real(8),parameter:: lightspeed=3.0d8 
    Real(8),parameter :: kB=1.3807d-23
    Real(8),parameter :: exp=2.71828d0
    Real(8),parameter :: ElectronCharge=1.6022d-19
    !Real(8),parameter :: ElectronMass=9.1095d-31
    Real(8),Parameter :: Epsilon=8.8542d-12
    Real(8),Parameter :: JtoeV=1.6022d-19
    Real(8),Parameter :: eVtoK=11605.d0
    Real(8),Parameter :: mTorrtoPa=0.13332d0
    Real(8),Parameter :: MinReal=1.d-15
    Real(8),save  :: R
    contains
    
             ! see example of RANDOM_SEED
    
           Subroutine DRandomInt(RR)
               ! USE IFPORT
                Implicit none
                 Real(8) :: RR
                 !Integer :: AA
                 !CALL init_random_seed()
                 Return 
           Entry  DRandom(RR)
                   CALL RANDOM_NUMBER(RR)
                   !Write(*,*) RR
                Return
          End  subroutine DRandomInt 
End  Module Constants

Module TypeModule
    Use ISO_C_BINDING
    Implicit none
    !  This Section is the definitions of the physicla parameters.
                Type  OneSpecy      ! Define a physical particle.
                     Integer(4) :: Model 
                     Real(8) :: Mass
                     Real(8) :: Vfactor
                End  Type  OneSpecy
                
                   !  ReactionType is definition of the reactions:
                   !   -1-The injected particle are removed because the generated particle is no longer included in simulations . 
                   !   0-no collision
                   !   1~99 for electrons  (1~10 Isotropic : 1-Elastic; 2-Exitation; 3- Ionization; 4-Attachment; 5-Dissociation)
                   !                                    (11-20 Anisotropic Ar: 11-Elastic; 12-Exitation; 13- Ionization;)
                   !   100~199 for  ions  (101~110 non-reactive : 101-isotropic; 102-Anisotropic; 103-ChargeExchange)
                   !                                    (111~120 reactive : 111-original particle; 112-new particle)    
                Type  ReactionOne    ! Define a collosion process.
                     Integer(4) :: Reactant,ReactionType,Resultant
                     Real(8) ::  Threshold
                End Type  ReactionOne
    
    
                !                struct MC_struct{
                !  double xp; // particle position                                                                               
                !  double yp;
                !  double zp;
                !
                !  double up; // particle velocity                                                                               
                !  double vp;
                !  double wp;
                !
                !  double qp; // the q                                                                                           
                !  unsigned long IDp;
                !
                !  int sp; // species of the particle                                                                            
                !}
                
    !  This  Section is the definitions of the  Particles used in the simulations 
                    !Integer(4),parameter :: NParMax=600000_4
                    !Integer(4),parameter :: NParMaxSmall=Int(0.9*dble(NParMax))
                !   One  Particle, one should change this to make 1D/2D and or Implicit/Explict simulations.  
                   Type, BIND (C) :: ParticleOne
                        Real(C_DOUBLE) ::  X
                        Real(C_DOUBLE) ::  Y
                        Real(C_DOUBLE) ::  Z
                        
                        Real(C_DOUBLE) ::  Vx
                        Real(C_DOUBLE) ::  Vy
                        Real(C_DOUBLE) ::  Vz
                        
                        Real(C_DOUBLE) ::  qp
                        Integer(C_LONG_LONG) ::  IDp
                        
                        
                        Integer(C_int) ::  sp_in
                        Integer(C_int) ::  sp_out
                        
                        Integer(C_int) :: CollisionType
                        
                        Real(C_DOUBLE) ::  dP_x
                        Real(C_DOUBLE) ::  dP_y
                        Real(C_DOUBLE) ::  dP_z
                        
                        Real(C_DOUBLE) ::  dE
                   EndType ParticleOne
              !   Particle  Bundles  for collision particle storage.
                  
                        !  This tpye defines the values for the collision praticle.
                         !  Index- ReactionType, same as the ReactionType defined below.
                         !  PhaseSpace- the position and velocity of the indection particle.
                         !  GasParticle- the position and velocity of the Gas particle .
                   Type MCCParticleOne
                          Integer(4) :: Index                            
                          Type(ParticleOne) :: PhaseSpace
                          Type(ParticleOne) :: Particle_Old
                          Real(8) :: Mass,V2,V,Energy,Energy_Old
                          Type(ParticleOne) ::  GasParticle
                    End Type  MCCParticleOne 
    
     ! This  Sections Defines the gas and collision properties used For PIC simulations. 
                       Type Gas
                             Integer(4) ::  NSpecy
                             Real(8) :: MGas=1.6726219d-27,TGas,NGas
                             Real(8) :: MFactor,EFactor
                       End Type  Gas
                           
                       Integer(4),parameter,private ::  NReactionMax=3_4
                       Integer(4),parameter,private :: NSigmaMax=150000_4 
                         
                       Type  ReactionBundle        !  Reactive bundle for all collision process of one particle withe one gas.
                             Integer(4) ::  NReaction 
                             Type(ReactionOne) ::  Reaction(1:NReactionMax)
                        End Type ReactionBundle
                       
                       !  Mcc bundle for all collision process of one particle withe one gas. 
                       ! Model is the index for the different MCC models:
                       ! 0-no collision; 1-elelctrons,the gas moculars are stationary; 2-nonreactive ions; 3-reactive ions;
                         
                        Type MCCBundle
                            Integer(4) ::  Model,NReaction=0,NSigma=0
                            Real(8) :: Mass,VFactor,CollisionRatio,EnergyInterval
                            Real(8) :: MFactor,EFactor
                            Real(8) :: Probility(1:NSigmaMax)
                            Type(ReactionOne) :: Reaction(0:NReactionMax)
                        End Type MCCBundle
                        
                        Integer(4),parameter,private :: NSigmaRaw=100_4
                              Type  SigmaRaw
                                   Real(8)  :: Threshold,MaxEnergy
                                   Integer(4) :: NSigma
                                   Real(8) :: Energy(NSigmaRaw),Sigma(NSigmaRaw)
                              end Type  SigmaRaw
                              
                              Type SigmaNormalized
                                  Real(8) :: EnergyInterval=0.d0
                                  Integer(4) :: NReaction=0,NSigma=0
                                  Real(8) ::  Value(1:NSigmaMax)
                              end Type SigmaNormalized
End Module TypeModule