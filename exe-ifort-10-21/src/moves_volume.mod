	  �/  �   k820309              13.0        "4wa                                                                                                           
       /home/rs/software/MCCCS-MN-10-21/src/moves_volume.F90 MOVES_VOLUME              VOLUME_1BOX VOLUME_2BOX INIT_MOVES_VOLUME UPDATE_VOLUME_MAX_DISPLACEMENT OUTPUT_VOLUME_STATS READ_CHECKPOINT_VOLUME WRITE_CHECKPOINT_VOLUME ALLOW_CUTOFF_FAILURE RESTORE_DISPL_TRANSL ACSVOL ACNVOL ACSHMAT ACNHMAT BSVOL BNVOL BSHMAT BNHMAT ACC_DISPL                      @                              
       RANDOM                      @                              
       ERR_EXIT                      @                              
       CONSTRUCT_KDTREE SCALE_KDTREE                                                    
                            @                              
                            @                              
       RECIP CALP SAVE_KVECTOR RESTORE_KVECTOR                      @                              
       SUMUP                @  !@                              '                    #VAL 	                �                              	                
   %         @     !                           
                  
       #UTIL_RANDOM!RANDOM%MT_STATE1    #UTIL_RANDOM!RANDOM%MT_STATE2    #UTIL_RANDOM!RANDOM%MT_MASK3    #UTIL_RANDOM!RANDOM%MT_MAG01    #RANDOM%ISHFT    #RANDOM%IOR    #RANDOM%IAND    #RANDOM%IEOR    #RANDOM%DBLE    #ISTREAM                   @                                                  #RANDOM%MT_STATE1%MTI    #RANDOM%MT_STATE1%INITIALIZED              �   @        �                                               �   @        �                                                   @                             �	                    #RANDOM%MT_STATE2%MT              �   @  �      �                      p                           p           & p         p o          p p                                               @                                                  #RANDOM%MT_MASK3%UPPER_MASK    #RANDOM%MT_MASK3%LOWER_MASK    #RANDOM%MT_MASK3%MATRIX_A    #RANDOM%MT_MASK3%T1_MASK    #RANDOM%MT_MASK3%T2_MASK              �   @        �                                               �   @        �                                              �   @        �                                              �   @        �                                              �   @        �                                                   @                                                  #RANDOM%MT_MAG01%MAG01              �   @  �      �                                                  p           & p         p            p                                                @                                 ISHFT               @                                 IOR               @                                 IAND               @                                 IEOR               @                                 DBLE           
   @                                         #         @       !                                              #ERR_EXIT%TRIM    #FILE     #LINENO !   #MSG "   #CODE #                 @                                 TRIM           
   @                                                  1           
   @                              !                     
   @                             "                    1           
   @                              #           #         @       !                           $                   #CONSTRUCT_KDTREE%REAL %   #CONSTRUCT_KDTREE%ASSOCIATED &   #CONSTRUCT_KDTREE%SIZE '   #CONSTRUCT_KDTREE%MOD (   #CONSTRUCT_KDTREE%FLOOR )   #CONSTRUCT_KDTREE%LOG *   #IBOX +   #ITREE ,   #LOUTPUT -                 @                            %     REAL               @                            &     ASSOCIATED               @                            '     SIZE               @                            (     MOD               @                            )     FLOOR               @                            *     LOG           
   @                              +                     
   @                              ,                       @                              -            #         @       !                           .                    #IBOX /   #FAC 0             
   @                              /                     
   @                              0     
      #         @       !                           1                   #RECIP%REAL 2   #RECIP%ABS 3   #RECIP%COS 4   #RECIP%SIN 5   #RECIP%MOD 6   #RECIP%INT 7   #IBOX 8   #VRECIPNEW 9   #VRECIPOLD :   #TYPE ;                 @              @             2     REAL               @                            3     ABS               @                            4     COS               @                            5     SIN               @                            6     MOD               @                            7     INT             @                              8                        @                              9     
                   @                              :     
                   @                              ;                     @    !                            <                   
                &                                           #         @       !                           =                    #IBOX >             
   @                              >           #         @       !                           ?                    #IBOX @             
   @                              @           #         @       !                           A                   #SUMUP%SQRT B   #SUMUP%ANY C   #SUMUP%INT D   #SUMUP%EXP E   #OVRLAP F   #V G   #IBOX H   #LVOL I                 @                            B     SQRT               @                            C     ANY               @                            D     INT               @                            E     EXP             @                              F                        @                              G                   
     p          p            p                                      @                              H                        @                              I                                                      J                   
                &                                                    @                                K                   
                &                   &                                           #         @                                  L                   #MATOPS%SQRT M   #MATOPS%ABS N   #MATOPS%ACOS O   #IBOX P                 @                            M     SQRT               @                            N     ABS               @                            O     ACOS           
   @                              P                                                     Q                   
                &                   &                                           $         @      !                         R                           #NUMBER S                     
   @                              S                    @                                T                   
                &                   &                                           #         @                                  U                    #BUILD_LINKED_CELL%INT V   #BUILD_LINKED_CELL%DBLE W   #BUILD_LINKED_CELL%ANINT X                 @                            V     INT               @                            W     DBLE               @                            X     ANINT #         @                                   Y                    #VOLUME_1BOX%LOG Z   #VOLUME_1BOX%EXP [   #VOLUME_1BOX%ABS \   #VOLUME_1BOX%MIN ]   #VOLUME_1BOX%MINVAL ^   #VOLUME_1BOX%SQRT _                 @                            Z     LOG               @                            [     EXP               @                            \     ABS               @                            ]     MIN               @                            ^     MINVAL               @                            _     SQRT #         @                                   `                    #VOLUME_2BOX%LOG a   #VOLUME_2BOX%EXP b   #VOLUME_2BOX%ABS c   #VOLUME_2BOX%MIN d   #VOLUME_2BOX%MINVAL e   #VOLUME_2BOX%SQRT f                 @                            a     LOG               @                            b     EXP               @                            c     ABS               @                            d     MIN               @                            e     MINVAL               @                            f     SQRT #         @                                   g                   #INIT_MOVES_VOLUME%ANY h   #INIT_MOVES_VOLUME%ALLOCATED i   #INIT_MOVES_VOLUME%ALL j   #INIT_MOVES_VOLUME%REAL k   #IO_INPUT l   #LPRINT m                 @                            h     ANY               @                            i     ALLOCATED               @                            j     ALL               @             @              k     REAL           
   @                              l                     
   @                              m           #         @                                   n                    #IO_OUTPUT o             
   @                              o           #         @                                   p                    #IO_OUTPUT q             
   @                              q           #         @                                   r                    #IO_CHKPT s             
   @                              s           #         @                                   t                    #IO_CHKPT u             
   @                              u                     D@@                               v            #         @                                  w                    #BOX x             
   @                              x                    @@                               y                   
                &                                                    @@                               z                   
                &                                                    @@                               {                   
                &                   &                                                    @@                               |                   
                &                   &                                                    @@                               }                   
                &                                                    @@                               ~                   
                &                                                    @@                                                  
                &                   &                                                    @@                               �                   
                &                   &                                                    @                                �                   
                &                                              �   K      fn#fn "   �     b   uapp(MOVES_VOLUME    �  G   J  UTIL_RANDOM    :  I   J  UTIL_RUNTIME    �  ^   J  UTIL_KDTREE    �  @   j  SIM_SYSTEM    !  @   J  SIM_CELL    a  h   J  ENERGY_KSPACE     �  F   J  ENERGY_PAIRWISE !     Y      REALPTR+VAR_TYPE %   h  H   a   REALPTR%VAL+VAR_TYPE #   �  8      RANDOM+UTIL_RANDOM C   �  �   �  UTIL_RANDOM!RANDOM%MT_STATE1+UTIL_RANDOM=MT_STATE1 5   t  H     RANDOM%MT_STATE1%MTI+UTIL_RANDOM=MTI E   �  H     RANDOM%MT_STATE1%INITIALIZED+UTIL_RANDOM=INITIALIZED C     i   �  UTIL_RANDOM!RANDOM%MT_STATE2+UTIL_RANDOM=MT_STATE2 3   m  �     RANDOM%MT_STATE2%MT+UTIL_RANDOM=MT A   !  �   �  UTIL_RANDOM!RANDOM%MT_MASK3+UTIL_RANDOM=MT_MASK3 B   		  H     RANDOM%MT_MASK3%UPPER_MASK+UTIL_RANDOM=UPPER_MASK B   Q	  H     RANDOM%MT_MASK3%LOWER_MASK+UTIL_RANDOM=LOWER_MASK >   �	  H     RANDOM%MT_MASK3%MATRIX_A+UTIL_RANDOM=MATRIX_A <   �	  H     RANDOM%MT_MASK3%T1_MASK+UTIL_RANDOM=T1_MASK <   )
  H     RANDOM%MT_MASK3%T2_MASK+UTIL_RANDOM=T2_MASK A   q
  k   �  UTIL_RANDOM!RANDOM%MT_MAG01+UTIL_RANDOM=MT_MAG01 8   �
  �     RANDOM%MT_MAG01%MAG01+UTIL_RANDOM=MAG01 /   �  >      RANDOM%ISHFT+UTIL_RANDOM=ISHFT +   �  <      RANDOM%IOR+UTIL_RANDOM=IOR -   
  =      RANDOM%IAND+UTIL_RANDOM=IAND -   G  =      RANDOM%IEOR+UTIL_RANDOM=IEOR -   �  =      RANDOM%DBLE+UTIL_RANDOM=DBLE +   �  @   e   RANDOM%ISTREAM+UTIL_RANDOM &     �       ERR_EXIT+UTIL_RUNTIME 0   �  =      ERR_EXIT%TRIM+UTIL_RUNTIME=TRIM +   �  L   e   ERR_EXIT%FILE+UTIL_RUNTIME -     @   e   ERR_EXIT%LINENO+UTIL_RUNTIME *   N  L   e   ERR_EXIT%MSG+UTIL_RUNTIME +   �  @   e   ERR_EXIT%CODE+UTIL_RUNTIME -   �        CONSTRUCT_KDTREE+UTIL_KDTREE 7   �  =      CONSTRUCT_KDTREE%REAL+UTIL_KDTREE=REAL C   (  C      CONSTRUCT_KDTREE%ASSOCIATED+UTIL_KDTREE=ASSOCIATED 7   k  =      CONSTRUCT_KDTREE%SIZE+UTIL_KDTREE=SIZE 5   �  <      CONSTRUCT_KDTREE%MOD+UTIL_KDTREE=MOD 9   �  >      CONSTRUCT_KDTREE%FLOOR+UTIL_KDTREE=FLOOR 5   "  <      CONSTRUCT_KDTREE%LOG+UTIL_KDTREE=LOG 2   ^  @   e   CONSTRUCT_KDTREE%IBOX+UTIL_KDTREE 3   �  @   e   CONSTRUCT_KDTREE%ITREE+UTIL_KDTREE 5   �  @   e   CONSTRUCT_KDTREE%LOUTPUT+UTIL_KDTREE )     [       SCALE_KDTREE+UTIL_KDTREE .   y  @   e   SCALE_KDTREE%IBOX+UTIL_KDTREE -   �  @   e   SCALE_KDTREE%FAC+UTIL_KDTREE $   �  �       RECIP+ENERGY_KSPACE .   �  =      RECIP%REAL+ENERGY_KSPACE=REAL ,     <      RECIP%ABS+ENERGY_KSPACE=ABS ,   G  <      RECIP%COS+ENERGY_KSPACE=COS ,   �  <      RECIP%SIN+ENERGY_KSPACE=SIN ,   �  <      RECIP%MOD+ENERGY_KSPACE=MOD ,   �  <      RECIP%INT+ENERGY_KSPACE=INT )   7  @   e   RECIP%IBOX+ENERGY_KSPACE .   w  @   e   RECIP%VRECIPNEW+ENERGY_KSPACE .   �  @   e   RECIP%VRECIPOLD+ENERGY_KSPACE )   �  @   e   RECIP%TYPE+ENERGY_KSPACE #   7  �       CALP+ENERGY_KSPACE +   �  R       SAVE_KVECTOR+ENERGY_KSPACE 0     @   e   SAVE_KVECTOR%IBOX+ENERGY_KSPACE .   U  R       RESTORE_KVECTOR+ENERGY_KSPACE 3   �  @   e   RESTORE_KVECTOR%IBOX+ENERGY_KSPACE &   �  �       SUMUP+ENERGY_PAIRWISE 0   �  =      SUMUP%SQRT+ENERGY_PAIRWISE=SQRT .   �  <      SUMUP%ANY+ENERGY_PAIRWISE=ANY .     <      SUMUP%INT+ENERGY_PAIRWISE=INT .   H  <      SUMUP%EXP+ENERGY_PAIRWISE=EXP -   �  @   e   SUMUP%OVRLAP+ENERGY_PAIRWISE (   �  �   e   SUMUP%V+ENERGY_PAIRWISE +   X  @   e   SUMUP%IBOX+ENERGY_PAIRWISE +   �  @   e   SUMUP%LVOL+ENERGY_PAIRWISE "   �  �       CELL_VOL+SIM_CELL    d  �       HMAT+SIM_CELL       �       MATOPS+SIM_CELL *   �  =      MATOPS%SQRT+SIM_CELL=SQRT (   �  <      MATOPS%ABS+SIM_CELL=ABS *     =      MATOPS%ACOS+SIM_CELL=ACOS %   B  @   e   MATOPS%IBOX+SIM_CELL #   �  �       MIN_WIDTH+SIM_CELL .   &  d       INTEGER_TO_STRING+UTIL_STRING 5   �  @   e   INTEGER_TO_STRING%NUMBER+UTIL_STRING    �  �       HMATI+SIM_CELL +   n  �       BUILD_LINKED_CELL+SIM_CELL 3   
   <      BUILD_LINKED_CELL%INT+SIM_CELL=INT 5   F   =      BUILD_LINKED_CELL%DBLE+SIM_CELL=DBLE 7   �   >      BUILD_LINKED_CELL%ANINT+SIM_CELL=ANINT    �   �       VOLUME_1BOX     �!  <      VOLUME_1BOX%LOG     �!  <      VOLUME_1BOX%EXP     "  <      VOLUME_1BOX%ABS     ?"  <      VOLUME_1BOX%MIN #   {"  ?      VOLUME_1BOX%MINVAL !   �"  =      VOLUME_1BOX%SQRT    �"  �       VOLUME_2BOX     �#  <      VOLUME_2BOX%LOG     �#  <      VOLUME_2BOX%EXP     9$  <      VOLUME_2BOX%ABS     u$  <      VOLUME_2BOX%MIN #   �$  ?      VOLUME_2BOX%MINVAL !   �$  =      VOLUME_2BOX%SQRT "   -%  �       INIT_MOVES_VOLUME &   &  <      INIT_MOVES_VOLUME%ANY ,   >&  B      INIT_MOVES_VOLUME%ALLOCATED &   �&  <      INIT_MOVES_VOLUME%ALL '   �&  =      INIT_MOVES_VOLUME%REAL +   �&  @   a   INIT_MOVES_VOLUME%IO_INPUT )   9'  @   a   INIT_MOVES_VOLUME%LPRINT /   y'  W       UPDATE_VOLUME_MAX_DISPLACEMENT 9   �'  @   a   UPDATE_VOLUME_MAX_DISPLACEMENT%IO_OUTPUT $   (  W       OUTPUT_VOLUME_STATS .   g(  @   a   OUTPUT_VOLUME_STATS%IO_OUTPUT '   �(  V       READ_CHECKPOINT_VOLUME 0   �(  @   a   READ_CHECKPOINT_VOLUME%IO_CHKPT (   =)  V       WRITE_CHECKPOINT_VOLUME 1   �)  @   a   WRITE_CHECKPOINT_VOLUME%IO_CHKPT %   �)  @       ALLOW_CUTOFF_FAILURE %   *  Q       RESTORE_DISPL_TRANSL )   d*  @   a   RESTORE_DISPL_TRANSL%BOX    �*  �       ACSVOL    0+  �       ACNVOL    �+  �       ACSHMAT    `,  �       ACNHMAT    -  �       BSVOL    �-  �       BNVOL    .  �       BSHMAT    �.  �       BNHMAT    d/  �       ACC_DISPL 