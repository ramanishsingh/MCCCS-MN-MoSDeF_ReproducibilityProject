	  �>  �   k820309              13.0        #4wa                                                                                                           
       /home/rs/software/MCCCS-MN-10-21/src/transfer_swap.F90 TRANSFER_SWAP       
       SWAP INIT_SWAP CNT OUTPUT_SWAP_STATS READ_CHECKPOINT_SWAP WRITE_CHECKPOINT_SWAP COMPUTE_BEG BEG ACCHEM BNCHEM                      @                              
       RANDOM                      @                              
       ERR_EXIT                                                    
                            @                              
                            @                              
       RECIP                      @                              
       ENERGY BOLTZ CORU                      @                              
       U_BONDED                      @                              
       ROSENBLUTH SCHEDULE EXPLCT SAFESCHEDULE                @  !@                         	     '                    #VAL 
                �                              
                
   #         @       !                                              #SCHEDULE%REAL    #SCHEDULE%INT    #IGROW    #IMOLTY    #INDEX    #IUTRY    #IPREV    #MOVETYPE    #GROUPTYPE                  @              @                  REAL               @                                 INT             @                                                      @                                                      @                                                      @                                                      @                                                      @                                                      @                                          %         @     !                                             
       #UTIL_RANDOM!RANDOM%MT_STATE1    #UTIL_RANDOM!RANDOM%MT_STATE2    #UTIL_RANDOM!RANDOM%MT_MASK3    #UTIL_RANDOM!RANDOM%MT_MAG01 !   #RANDOM%ISHFT #   #RANDOM%IOR $   #RANDOM%IAND %   #RANDOM%IEOR &   #RANDOM%DBLE '   #ISTREAM (                  @                                                  #RANDOM%MT_STATE1%MTI    #RANDOM%MT_STATE1%INITIALIZED              �   @        �                                               �   @        �                                                   @                             �	                    #RANDOM%MT_STATE2%MT              �   @  �      �                      p                           p           & p         p o          p p                                               @                                                  #RANDOM%MT_MASK3%UPPER_MASK    #RANDOM%MT_MASK3%LOWER_MASK    #RANDOM%MT_MASK3%MATRIX_A    #RANDOM%MT_MASK3%T1_MASK    #RANDOM%MT_MASK3%T2_MASK               �   @        �                                               �   @        �                                              �   @        �                                              �   @        �                                              �   @        �                                                    @                        !                          #RANDOM%MT_MAG01%MAG01 "             �   @  �      �                 "                                 p           & p         p            p                                                @                            #     ISHFT               @                            $     IOR               @                            %     IAND               @                            &     IEOR               @                            '     DBLE           
   @                              (           #         @       !                           )                   #ERR_EXIT%TRIM *   #FILE +   #LINENO ,   #MSG -   #CODE .                 @                            *     TRIM           
   @                             +                    1           
   @                              ,                     
   @                             -                    1           
   @                              .           #         @       !                           /                   #RECIP%REAL 0   #RECIP%ABS 1   #RECIP%COS 2   #RECIP%SIN 3   #RECIP%MOD 4   #RECIP%INT 5   #IBOX 6   #VRECIPNEW 7   #VRECIPOLD 8   #TYPE 9                 @              @             0     REAL               @                            1     ABS               @                            2     COS               @                            3     SIN               @                            4     MOD               @                            5     INT             @                              6                        @                              7     
                   @                              8     
                   @                              9            #         @       !                           :                   #ENERGY%SQRT ;   #ENERGY%INT <   #ENERGY%EXP =   #I >   #IMOLTY ?   #V @   #FLAGON A   #IBOX B   #ISTART C   #IUEND D   #LLJII E   #OVRLAP F   #LTORS G   #LCHARGE_TABLE H   #LFAVOR I   #LATOM_TRAXYZ J                 @                            ;     SQRT               @                            <     INT               @                            =     EXP             @                              >                        @                              ?                        @                              @                   
     p          p            p                                      @                              A                        @                              B                        @                              C                        @                              D                        @                              E                        @                              F                        @                              G                        @                              H                        @                              I                        @                              J            #         @       !                           K                  #NUMAX L   #NCHMAX M   #BOLTZ%SQRT N   #BOLTZ%ANY O   #BOLTZ%INT P   #BOLTZ%EXP Q   #LNEW R   #LFIRST S   #OVRLAP T   #I U   #ICHARGE V   #IMOLTY W   #IBOX X   #ICHOI Y   #IUFROM Z   #NTOGROW [   #GLIST \   #MAXLEN ]                @                             L                         @                             M                          @                            N     SQRT               @                            O     ANY               @                            P     INT               @                            Q     EXP             @                              R                        @                              S                        @                              T                        @                              U                        @                              V                        @                              W                        @                              X                        @                              Y                        @                              Z                        @                              [                       @                              \                         p          5 r L       5 r L                                 @                              ]     
       %         @     !                           ^                   
       #CORU%ALL _   #CORU%SQRT `   #CORU%EXP a   #IMOLTY b   #JMOLTY c   #RHO d   #IBOX e                 @                            _     ALL               @                            `     SQRT               @                            a     EXP           
   @                              b                     
   @                              c                     
   @                              d     
                
   @                              e           #         @       !                �          f                   #U_BONDED%SQRT g   #U_BONDED%ACOS h   #I i   #IMOLTY j   #VVIB k   #VBEND l   #VTG m                 @                            g     SQRT               @                            h     ACOS           
   @                              i                     
   @                              j                       @                              k     
                   @                              l     
                   @                              m     
       #         @       !                           n                   #ROSENBLUTH%REAL o   #ROSENBLUTH%EXP p   #ROSENBLUTH%ANINT q   #ROSENBLUTH%SQRT r   #ROSENBLUTH%SUM s   #LNEW t   #LTERM u   #I v   #ICHARGE w   #IMOLTY x   #IFROM y   #IBOX z   #IGROW {   #WADD |   #LFIXNOW }   #CWTORF ~   #MOVETYPE    #FIRST_BEAD �   #SECOND_BEAD �                 @              @             o     REAL               @                            p     EXP               @                            q     ANINT               @                            r     SQRT               @                            s     SUM             @                              t                        @                              u                        @                              v                        @                              w                        @                              x                        @                              y                        @                              z                        @                              {                        @                              |     
                   @                              }                        @                              ~     
                   @                                                      @                              �                        @                              �            #         @       !                           �                   #EXPLCT%EXP �   #EXPLCT%SQRT �   #EXPLCT%SIN �   #EXPLCT%COS �   #ICHAIN �   #VMETHYL �   #LCRYSL �   #LSWITCH �                 @                            �     EXP               @                            �     SQRT               @                            �     SIN               @                            �     COS             @                              �                        @                              �     
                   @                              �                        @                              �            #         @       !                           �                   #SAFESCHEDULE%INT �   #SAFESCHEDULE%DBLE �   #SAFESCHEDULE%PRESENT �   #IGROW �   #IMOLTY �   #ISLEN �   #IUTRY �   #FINDEX �   #MOVETYPE �   #IPREV �                 @                            �     INT               @                            �     DBLE               @                            �     PRESENT             @                              �                        @                              �                        @                              �                        @                              �                        @                              �                        @                              �                        @                              �            #         @                                 �                   #MIMAGE%AINT �   #MIMAGE%SIGN �   #RXUIJ �   #RYUIJ �   #RZUIJ �   #IBOX �                 @                            �     AINT               @                            �     SIGN           
  @                              �     
                 
  @                              �     
                 
  @                              �     
                 
   @                              �           #         @                                  �                    #IBOX �             
   @                              �                                                     �                   
                &                   &                                                                                     �                   
                &                                           #         @                                  �                   #UPDATE_LINKED_CELL%INT �   #IMOL �                 @                            �     INT           
   @                              �           #         @                                   �                    #SWAP%LOG �   #SWAP%AINT �   #SWAP%LOG10 �   #SWAP%EXP �   #SWAP%SQRT �   #SWAP%DBLE �   #SWAP%INT �   #SWAP%REAL �                 @                            �     LOG               @                            �     AINT               @                            �     LOG10               @                            �     EXP               @                            �     SQRT               @                            �     DBLE               @                            �     INT               @             @              �     REAL #         @                                   �                   #INIT_SWAP%ALLOCATED �   #INIT_SWAP%REAL �   #IO_INPUT �   #LPRINT �                                                      @                            �     ALLOCATED               @             @              �     REAL           
  @@                              �                     
   @                              �           #         @                                   �                    #CNT%DBLE �                 @                            �     DBLE #         @                                   �                    #IO_OUTPUT �             
   @                              �           #         @                                   �                    #IO_CHKPT �             
   @                              �           #         @                                   �                    #IO_CHKPT �             
   @                              �           #         @                                  �                   #COMPUTE_BEG%INT �   #COMPUTE_BEG%REAL �   #IMOLTY �                 @                            �     INT               @             @              �     REAL           
   @                              �                      @                                �                     @@                               �                   
                &                   &                                                    @@                               �                                   &                   &                                              �   M      fn#fn #   �   ~   b   uapp(TRANSFER_SWAP    k  G   J  UTIL_RANDOM    �  I   J  UTIL_RUNTIME    �  @   j  SIM_SYSTEM    ;  @   J  SIM_CELL    {  F   J  ENERGY_KSPACE     �  R   J  ENERGY_PAIRWISE &     I   J  ENERGY_INTRAMOLECULAR    \  h   J  MOVES_CBMC !   �  Y      REALPTR+VAR_TYPE %     H   a   REALPTR%VAL+VAR_TYPE $   e  �       SCHEDULE+MOVES_CBMC .   '  =      SCHEDULE%REAL+MOVES_CBMC=REAL ,   d  <      SCHEDULE%INT+MOVES_CBMC=INT *   �  @   e   SCHEDULE%IGROW+MOVES_CBMC +   �  @   e   SCHEDULE%IMOLTY+MOVES_CBMC *      @   e   SCHEDULE%INDEX+MOVES_CBMC *   `  @   e   SCHEDULE%IUTRY+MOVES_CBMC *   �  @   e   SCHEDULE%IPREV+MOVES_CBMC -   �  @   e   SCHEDULE%MOVETYPE+MOVES_CBMC .      @   e   SCHEDULE%GROUPTYPE+MOVES_CBMC #   `  8      RANDOM+UTIL_RANDOM C   �  �   �  UTIL_RANDOM!RANDOM%MT_STATE1+UTIL_RANDOM=MT_STATE1 5   $	  H     RANDOM%MT_STATE1%MTI+UTIL_RANDOM=MTI E   l	  H     RANDOM%MT_STATE1%INITIALIZED+UTIL_RANDOM=INITIALIZED C   �	  i   �  UTIL_RANDOM!RANDOM%MT_STATE2+UTIL_RANDOM=MT_STATE2 3   
  �     RANDOM%MT_STATE2%MT+UTIL_RANDOM=MT A   �
  �   �  UTIL_RANDOM!RANDOM%MT_MASK3+UTIL_RANDOM=MT_MASK3 B   �  H     RANDOM%MT_MASK3%UPPER_MASK+UTIL_RANDOM=UPPER_MASK B     H     RANDOM%MT_MASK3%LOWER_MASK+UTIL_RANDOM=LOWER_MASK >   I  H     RANDOM%MT_MASK3%MATRIX_A+UTIL_RANDOM=MATRIX_A <   �  H     RANDOM%MT_MASK3%T1_MASK+UTIL_RANDOM=T1_MASK <   �  H     RANDOM%MT_MASK3%T2_MASK+UTIL_RANDOM=T2_MASK A   !  k   �  UTIL_RANDOM!RANDOM%MT_MAG01+UTIL_RANDOM=MT_MAG01 8   �  �     RANDOM%MT_MAG01%MAG01+UTIL_RANDOM=MAG01 /   @  >      RANDOM%ISHFT+UTIL_RANDOM=ISHFT +   ~  <      RANDOM%IOR+UTIL_RANDOM=IOR -   �  =      RANDOM%IAND+UTIL_RANDOM=IAND -   �  =      RANDOM%IEOR+UTIL_RANDOM=IEOR -   4  =      RANDOM%DBLE+UTIL_RANDOM=DBLE +   q  @   e   RANDOM%ISTREAM+UTIL_RANDOM &   �  �       ERR_EXIT+UTIL_RUNTIME 0   5  =      ERR_EXIT%TRIM+UTIL_RUNTIME=TRIM +   r  L   e   ERR_EXIT%FILE+UTIL_RUNTIME -   �  @   e   ERR_EXIT%LINENO+UTIL_RUNTIME *   �  L   e   ERR_EXIT%MSG+UTIL_RUNTIME +   J  @   e   ERR_EXIT%CODE+UTIL_RUNTIME $   �  �       RECIP+ENERGY_KSPACE .   _  =      RECIP%REAL+ENERGY_KSPACE=REAL ,   �  <      RECIP%ABS+ENERGY_KSPACE=ABS ,   �  <      RECIP%COS+ENERGY_KSPACE=COS ,     <      RECIP%SIN+ENERGY_KSPACE=SIN ,   P  <      RECIP%MOD+ENERGY_KSPACE=MOD ,   �  <      RECIP%INT+ENERGY_KSPACE=INT )   �  @   e   RECIP%IBOX+ENERGY_KSPACE .     @   e   RECIP%VRECIPNEW+ENERGY_KSPACE .   H  @   e   RECIP%VRECIPOLD+ENERGY_KSPACE )   �  @   e   RECIP%TYPE+ENERGY_KSPACE '   �        ENERGY+ENERGY_PAIRWISE 1   �  =      ENERGY%SQRT+ENERGY_PAIRWISE=SQRT /     <      ENERGY%INT+ENERGY_PAIRWISE=INT /   T  <      ENERGY%EXP+ENERGY_PAIRWISE=EXP )   �  @   e   ENERGY%I+ENERGY_PAIRWISE .   �  @   e   ENERGY%IMOLTY+ENERGY_PAIRWISE )     �   e   ENERGY%V+ENERGY_PAIRWISE .   �  @   e   ENERGY%FLAGON+ENERGY_PAIRWISE ,   �  @   e   ENERGY%IBOX+ENERGY_PAIRWISE .   $  @   e   ENERGY%ISTART+ENERGY_PAIRWISE -   d  @   e   ENERGY%IUEND+ENERGY_PAIRWISE -   �  @   e   ENERGY%LLJII+ENERGY_PAIRWISE .   �  @   e   ENERGY%OVRLAP+ENERGY_PAIRWISE -   $  @   e   ENERGY%LTORS+ENERGY_PAIRWISE 5   d  @   e   ENERGY%LCHARGE_TABLE+ENERGY_PAIRWISE .   �  @   e   ENERGY%LFAVOR+ENERGY_PAIRWISE 4   �  @   e   ENERGY%LATOM_TRAXYZ+ENERGY_PAIRWISE &   $  #      BOLTZ+ENERGY_PAIRWISE !   G  @     NUMAX+SIM_SYSTEM "   �  @     NCHMAX+SIM_SYSTEM 0   �  =      BOLTZ%SQRT+ENERGY_PAIRWISE=SQRT .     <      BOLTZ%ANY+ENERGY_PAIRWISE=ANY .   @  <      BOLTZ%INT+ENERGY_PAIRWISE=INT .   |  <      BOLTZ%EXP+ENERGY_PAIRWISE=EXP +   �  @   e   BOLTZ%LNEW+ENERGY_PAIRWISE -   �  @   e   BOLTZ%LFIRST+ENERGY_PAIRWISE -   8  @   e   BOLTZ%OVRLAP+ENERGY_PAIRWISE (   x  @   e   BOLTZ%I+ENERGY_PAIRWISE .   �  @   e   BOLTZ%ICHARGE+ENERGY_PAIRWISE -   �  @   e   BOLTZ%IMOLTY+ENERGY_PAIRWISE +   8  @   e   BOLTZ%IBOX+ENERGY_PAIRWISE ,   x  @   e   BOLTZ%ICHOI+ENERGY_PAIRWISE -   �  @   e   BOLTZ%IUFROM+ENERGY_PAIRWISE .   �  @   e   BOLTZ%NTOGROW+ENERGY_PAIRWISE ,   8  �   e   BOLTZ%GLIST+ENERGY_PAIRWISE -   �  @   e   BOLTZ%MAXLEN+ENERGY_PAIRWISE %      �       CORU+ENERGY_PAIRWISE -   �   <      CORU%ALL+ENERGY_PAIRWISE=ALL /   �   =      CORU%SQRT+ENERGY_PAIRWISE=SQRT -   +!  <      CORU%EXP+ENERGY_PAIRWISE=EXP ,   g!  @   e   CORU%IMOLTY+ENERGY_PAIRWISE ,   �!  @   e   CORU%JMOLTY+ENERGY_PAIRWISE )   �!  @   e   CORU%RHO+ENERGY_PAIRWISE *   '"  @   e   CORU%IBOX+ENERGY_PAIRWISE /   g"  �       U_BONDED+ENERGY_INTRAMOLECULAR 9   #  =      U_BONDED%SQRT+ENERGY_INTRAMOLECULAR=SQRT 9   C#  =      U_BONDED%ACOS+ENERGY_INTRAMOLECULAR=ACOS 1   �#  @   e   U_BONDED%I+ENERGY_INTRAMOLECULAR 6   �#  @   e   U_BONDED%IMOLTY+ENERGY_INTRAMOLECULAR 4    $  @   e   U_BONDED%VVIB+ENERGY_INTRAMOLECULAR 5   @$  @   e   U_BONDED%VBEND+ENERGY_INTRAMOLECULAR 3   �$  @   e   U_BONDED%VTG+ENERGY_INTRAMOLECULAR &   �$  W      ROSENBLUTH+MOVES_CBMC 0   &  =      ROSENBLUTH%REAL+MOVES_CBMC=REAL .   T&  <      ROSENBLUTH%EXP+MOVES_CBMC=EXP 2   �&  >      ROSENBLUTH%ANINT+MOVES_CBMC=ANINT 0   �&  =      ROSENBLUTH%SQRT+MOVES_CBMC=SQRT .   '  <      ROSENBLUTH%SUM+MOVES_CBMC=SUM +   G'  @   e   ROSENBLUTH%LNEW+MOVES_CBMC ,   �'  @   e   ROSENBLUTH%LTERM+MOVES_CBMC (   �'  @   e   ROSENBLUTH%I+MOVES_CBMC .   (  @   e   ROSENBLUTH%ICHARGE+MOVES_CBMC -   G(  @   e   ROSENBLUTH%IMOLTY+MOVES_CBMC ,   �(  @   e   ROSENBLUTH%IFROM+MOVES_CBMC +   �(  @   e   ROSENBLUTH%IBOX+MOVES_CBMC ,   )  @   e   ROSENBLUTH%IGROW+MOVES_CBMC +   G)  @   e   ROSENBLUTH%WADD+MOVES_CBMC .   �)  @   e   ROSENBLUTH%LFIXNOW+MOVES_CBMC -   �)  @   e   ROSENBLUTH%CWTORF+MOVES_CBMC /   *  @   e   ROSENBLUTH%MOVETYPE+MOVES_CBMC 1   G*  @   e   ROSENBLUTH%FIRST_BEAD+MOVES_CBMC 2   �*  @   e   ROSENBLUTH%SECOND_BEAD+MOVES_CBMC "   �*  �       EXPLCT+MOVES_CBMC *   �+  <      EXPLCT%EXP+MOVES_CBMC=EXP ,   �+  =      EXPLCT%SQRT+MOVES_CBMC=SQRT *   �+  <      EXPLCT%SIN+MOVES_CBMC=SIN *   7,  <      EXPLCT%COS+MOVES_CBMC=COS )   s,  @   e   EXPLCT%ICHAIN+MOVES_CBMC *   �,  @   e   EXPLCT%VMETHYL+MOVES_CBMC )   �,  @   e   EXPLCT%LCRYSL+MOVES_CBMC *   3-  @   e   EXPLCT%LSWITCH+MOVES_CBMC (   s-  �       SAFESCHEDULE+MOVES_CBMC 0   T.  <      SAFESCHEDULE%INT+MOVES_CBMC=INT 2   �.  =      SAFESCHEDULE%DBLE+MOVES_CBMC=DBLE 8   �.  @      SAFESCHEDULE%PRESENT+MOVES_CBMC=PRESENT .   /  @   e   SAFESCHEDULE%IGROW+MOVES_CBMC /   M/  @   e   SAFESCHEDULE%IMOLTY+MOVES_CBMC .   �/  @   e   SAFESCHEDULE%ISLEN+MOVES_CBMC .   �/  @   e   SAFESCHEDULE%IUTRY+MOVES_CBMC /   0  @   e   SAFESCHEDULE%FINDEX+MOVES_CBMC 1   M0  @   e   SAFESCHEDULE%MOVETYPE+MOVES_CBMC .   �0  @   e   SAFESCHEDULE%IPREV+MOVES_CBMC     �0  �       MIMAGE+SIM_CELL *   b1  =      MIMAGE%AINT+SIM_CELL=AINT *   �1  =      MIMAGE%SIGN+SIM_CELL=SIGN &   �1  @   e   MIMAGE%RXUIJ+SIM_CELL &   2  @   e   MIMAGE%RYUIJ+SIM_CELL &   \2  @   e   MIMAGE%RZUIJ+SIM_CELL %   �2  @   e   MIMAGE%IBOX+SIM_CELL     �2  R       SETPBC+SIM_CELL %   .3  @   e   SETPBC%IBOX+SIM_CELL    n3  �       HMAT+SIM_CELL "   4  �       CELL_VOL+SIM_CELL ,   �4  n       UPDATE_LINKED_CELL+SIM_CELL 4   5  <      UPDATE_LINKED_CELL%INT+SIM_CELL=INT 1   H5  @   e   UPDATE_LINKED_CELL%IMOL+SIM_CELL    �5  �       SWAP    F6  <      SWAP%LOG    �6  =      SWAP%AINT    �6  >      SWAP%LOG10    �6  <      SWAP%EXP    97  =      SWAP%SQRT    v7  =      SWAP%DBLE    �7  <      SWAP%INT    �7  =      SWAP%REAL    ,8  �       INIT_SWAP $   �8  B      INIT_SWAP%ALLOCATED    "9  =      INIT_SWAP%REAL #   _9  @   a   INIT_SWAP%IO_INPUT !   �9  @   a   INIT_SWAP%LPRINT    �9  V       CNT    5:  =      CNT%DBLE "   r:  W       OUTPUT_SWAP_STATS ,   �:  @   a   OUTPUT_SWAP_STATS%IO_OUTPUT %   	;  V       READ_CHECKPOINT_SWAP .   _;  @   a   READ_CHECKPOINT_SWAP%IO_CHKPT &   �;  V       WRITE_CHECKPOINT_SWAP /   �;  @   a   WRITE_CHECKPOINT_SWAP%IO_CHKPT    5<         COMPUTE_BEG     �<  <      COMPUTE_BEG%INT !   �<  =      COMPUTE_BEG%REAL #   -=  @   a   COMPUTE_BEG%IMOLTY    m=  @       BEG    �=  �       ACCHEM    Q>  �       BNCHEM 