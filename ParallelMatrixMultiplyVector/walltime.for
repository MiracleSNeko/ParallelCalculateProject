C***********************************************************************
C
C     walltime.for
C     计算串行程序和并行程序的运行时间并显�?
C
C***********************************************************************

      PROGRAM WALLTIME_FOR
      
      REAL*4 T, TS, T1, T4(4), T8(8), T16(16)
            
      OPEN(8, FILE = 'walltime', ACCESS = 'DIRECT', FORM = 'UNFORMATTED
     & ', RECL = 4)
      READ(8, REC = 1) T
      CLOSE(8)
      PRINT *, "serial program walltime: ", T, "ms"
      
      OPEN(8, FILE = 'walltime_mpiio', ACCESS = 'DIRECT', 
     & FORM = 'UNFORMATTED', RECL = 4)
      READ(8, REC = 1) TS
      CLOSE(8)
      PRINT *, "serial program(MPI I/O) walltime: ", TS, "ms"
      
      OPEN(8, FILE = 'walltime1', ACCESS = 'DIRECT', FORM = 'UNFORMATTED
     & ', RECL = 4)
      READ(8, REC = 1) T1
      CLOSE(8)
      PRINT *, "parallel program walltime (1 process): ", T1, "ms"
      
      OPEN(8, FILE = 'walltime4', ACCESS = 'DIRECT', FORM = 'UNFORMATTED
     & ', RECL = 4*4)
      READ(8, REC = 1) T4
      CLOSE(8)
      SUM4 = .0
      DO 40 I = 1, 4
        SUM4 = SUM4 + T4(I)
   40 CONTINUE
      PRINT *, "parallel program walltime (4 process): ", SUM4/4, "ms"
      
      OPEN(8, FILE = 'walltime8', ACCESS = 'DIRECT', FORM = 'UNFORMATTED
     & ', RECL = 4*8)
      READ(8, REC = 1) T8
      CLOSE(8)
      SUM8 = .0
      DO 80 I = 1, 8
        SUM8 = SUM8 + T8(I)
   80 CONTINUE
      PRINT *, "parallel program walltime (8 process): ", SUM8/8, "ms"
      
      OPEN(8, FILE = 'walltime16', ACCESS = 'DIRECT', FORM = 
     & 'UNFORMATTED', RECL = 4*16)
      READ(8, REC = 1) T16
      CLOSE(8)
      SUM16 = .0
      DO 160 I = 1, 16
        SUM16 = SUM16 + T16(I)
  160 CONTINUE
      PRINT *, "parallel program walltime (16 process): ", SUM16/16, 
     & "ms"
      
      END PROGRAM WALLTIME_FOR