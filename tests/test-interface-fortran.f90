program test
   use gtoint
   implicit none

   real(8), parameter :: D_BOHR_A = 1.889726124559

   integer, parameter :: N0 = 4, N1 = 5, ND = 2

   real(8), parameter :: REF_OI(N0,N1,ND) = reshape((/ &
   & -0.235774,   0.129365,  -0.129365,   0.000000, &
   &  0.000000,   0.000000,   0.000000,  -0.205919, &
   & -0.136124,   0.074689,  -0.074689,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.205919, &
   &  0.000000,   0.205919,   0.205919,   0.000000, &
   &  0.041536,   0.437539,   0.135246,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.354850, &
   &  0.023981,   0.126788,  -0.047741,   0.000000, &
   &  0.000000,   0.000000,   0.000000,  -0.136914, &
   &  0.249533,  -0.491764,   0.000000,   0.000000  &
   & /), (/N0,N1,ND/))
   real(8), parameter :: REF_KEI(N0,N1,ND) = reshape((/ &
   & -0.866968,   1.029690,  -1.029690,   0.000000, &
   &  0.000000,   0.000000,   0.000000,  -0.650850, &
   & -0.500544,   0.594492,  -0.594492,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.650850, &
   &  0.000000,   0.650850,   0.650850,   0.000000, &
   &  0.756228,   0.792574,   1.674862,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   1.778606, &
   &  0.436608,   0.059897,   0.569286,   0.000000, &
   &  0.000000,   0.000000,   0.000000,  -1.089777, &
   &  0.917560,  -2.868384,  -0.000000,   0.000000  &
   & /), (/N0,N1,ND/))
   real(8), parameter :: REF_NAI(N0,N1,ND) = reshape((/ &
   & -0.401990,   0.232452,  -0.232452,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,  -0.253116, &
   & -0.232089,   0.134206,  -0.134206,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,   0.253116, &
   &  0.000000,   0.253116,   0.253116,   0.000000, &
   &  0.043037,   0.941273,   0.208513,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,   0.370920, &
   &  0.024847,   0.426718,   0.003659,  -0.000000, &
   &  0.000000,   0.000000,   0.000000,  -0.168744, &
   &  0.224916,  -0.539665,  -0.000000,  -0.000000  &
   & /), (/N0,N1,ND/))
   real(8), parameter :: REF_ECP(N0,N1,ND) = reshape((/ &
   &  0.000000,   0.018598,  -0.018598,   0.000000, &
   &  0.000000,   0.000000,   0.000000,  -0.009753, &
   &  0.000000,   0.010738,  -0.010738,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.009753, &
   &  0.000000,   0.009753,   0.009753,   0.000000, &
   & -0.118920,   0.000000,   0.000000,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.000000, &
   & -0.068658,   0.000000,   0.000000,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.000000, &
   & -0.061982,   0.000000,   0.000000,   0.000000  &
   & /), (/N0,N1,ND/))
   real(8), parameter :: REF_ERI(N0,N1,N1,N0,ND) = reshape((/ &
   &  0.107024,  -0.076361,   0.076361,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.069210, &
   &  0.058498,  -0.040842,   0.040842,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.069210, &
   &  0.000000,  -0.065501,  -0.065501,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.018863, &
   &  0.022965,  -0.025076,   0.012967,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.001191, &
   & -0.017261,   0.011056,  -0.015746,  -0.000000, &
   &  0.000000,  -0.000000,  -0.000000,   0.003709, &
   &  0.058498,  -0.040842,   0.040842,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.035675, &
   &  0.039477,  -0.029201,   0.029201,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.035675, &
   &  0.000000,  -0.042100,  -0.042100,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.018863, &
   & -0.017261,   0.015746,  -0.011056,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.001191, &
   &  0.022965,  -0.012967,   0.025076,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.003709, &
   &  0.000000,  -0.010463,  -0.010463,   0.000000, &
   &  0.000000,  -0.000000,  -0.000000,   0.008400, &
   &  0.000000,  -0.015740,  -0.015740,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.008400, &
   &  0.040226,  -0.032422,   0.032422,   0.000000, &
   & -0.076361,   0.067455,  -0.054466,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.047595, &
   & -0.040842,   0.036871,  -0.026752,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.044206, &
   & -0.010463,   0.049659,   0.035811,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.016858, &
   & -0.025076,   0.029815,  -0.014004,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.000218, &
   &  0.015746,  -0.011429,   0.014337,   0.000000, &
   & -0.000000,   0.000000,   0.000000,  -0.006451, &
   & -0.040842,   0.036871,  -0.026752,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.025605, &
   & -0.029201,   0.031454,  -0.017000,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.020086, &
   & -0.015740,   0.040060,   0.016597,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.010385, &
   &  0.011056,  -0.011429,   0.007273,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.000043, &
   & -0.012967,   0.010089,  -0.014004,  -0.000000, &
   &  0.000000,  -0.000000,  -0.000000,   0.000121, &
   & -0.065501,   0.049659,  -0.035811,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.054497, &
   & -0.042100,   0.040060,  -0.016597,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.041102, &
   & -0.032422,   0.082466,   0.026277,  -0.000000, &
   &  0.076361,  -0.054466,   0.067455,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.044206, &
   &  0.040842,  -0.026752,   0.036871,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.047595, &
   & -0.010463,  -0.035811,  -0.049659,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.010385, &
   &  0.012967,  -0.014004,   0.010089,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.000043, &
   & -0.011056,   0.007273,  -0.011429,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.000121, &
   &  0.040842,  -0.026752,   0.036871,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.020086, &
   &  0.029201,  -0.017000,   0.031454,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.025605, &
   & -0.015740,  -0.016597,  -0.040060,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.016858, &
   & -0.015746,   0.014337,  -0.011429,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.000218, &
   &  0.025076,  -0.014004,   0.029815,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.006451, &
   & -0.065501,   0.035811,  -0.049659,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.041102, &
   & -0.042100,   0.016597,  -0.040060,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.054497, &
   &  0.032422,   0.026277,   0.082466,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.022851, &
   &  0.018863,  -0.016858,   0.010385,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.004425, &
   & -0.018863,   0.010385,  -0.016858,  -0.000000, &
   &  0.000000,  -0.000000,  -0.000000,  -0.000000, &
   &  0.069210,  -0.047595,   0.044206,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.057033, &
   &  0.035675,  -0.025605,   0.020086,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.051709, &
   &  0.008400,  -0.054497,  -0.041102,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.004425, &
   &  0.001191,  -0.000218,   0.000043,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.004592, &
   & -0.001191,   0.000043,  -0.000218,  -0.000000, &
   &  0.000000,  -0.000000,  -0.000000,  -0.000000, &
   & -0.069210,   0.044206,  -0.047595,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.051709, &
   & -0.035675,   0.020086,  -0.025605,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.057033, &
   &  0.008400,   0.041102,   0.054497,   0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.000000, &
   &  0.003709,  -0.006451,  -0.000121,   0.000000, &
   & -0.000000,  -0.000000,  -0.000000,   0.000000, &
   &  0.003709,   0.000121,   0.006451,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.006572, &
   & -0.039609,  -0.155191,  -0.089850,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.114152, &
   & -0.017304,  -0.061048,  -0.015131,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.058507, &
   & -0.069739,   0.161476,  -0.000000,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.021638, &
   & -0.036174,   0.017621,  -0.026735,  -0.000000, &
   & -0.000000,  -0.000000,   0.000000,   0.003601, &
   &  0.013564,   0.016446,   0.018737,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.019417, &
   & -0.017304,  -0.091358,  -0.045441,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.059449, &
   & -0.019628,  -0.028005,  -0.015685,  -0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.027322, &
   & -0.046379,   0.106141,  -0.000000,  -0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.018491, &
   &  0.021241,  -0.000246,   0.018737,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.005418, &
   & -0.017906,  -0.021438,  -0.026735,  -0.000000, &
   &  0.000000,  -0.000000,  -0.000000,   0.008234, &
   & -0.000957,   0.012254,  -0.000000,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.027763, &
   & -0.015532,   0.039263,  -0.000000,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.000112, &
   & -0.044442,   0.017642,  -0.017821,  -0.000000, &
   &  0.036494,   0.084770,   0.077988,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.085171, &
   &  0.016427,   0.030462,   0.019251,   0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.046004, &
   &  0.052692,  -0.109091,   0.000121,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.022880, &
   &  0.043999,  -0.028961,   0.031686,   0.000000, &
   &  0.000000,   0.000000,  -0.000000,  -0.003522, &
   & -0.013821,  -0.013246,  -0.019258,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.023592, &
   &  0.014015,   0.053176,   0.038532,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.048925, &
   &  0.024779,  -0.001771,   0.019763,   0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.020110, &
   &  0.047043,  -0.081392,   0.002442,   0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.014710, &
   & -0.015078,   0.002878,  -0.014123,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.002038, &
   &  0.013953,   0.006626,   0.020291,   0.000000, &
   & -0.000000,   0.000000,   0.000000,  -0.005369, &
   &  0.018706,   0.096123,   0.044622,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.105566, &
   &  0.028873,  -0.001415,   0.002319,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.035093, &
   &  0.102579,  -0.155203,   0.022279,   0.000000, &
   & -0.033064,  -0.105916,  -0.077988,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.072985, &
   & -0.011967,  -0.047009,  -0.019251,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.041460, &
   & -0.039134,   0.115809,   0.000121,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.012354, &
   & -0.020686,   0.009529,  -0.020291,  -0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.003398, &
   &  0.009375,   0.009675,   0.014123,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.007169, &
   & -0.014379,  -0.062859,  -0.038532,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.030424, &
   & -0.011994,  -0.031007,  -0.019763,  -0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.018424, &
   & -0.016804,   0.072380,   0.002442,  -0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.016624, &
   &  0.019846,  -0.000699,   0.019258,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.007134, &
   & -0.020243,  -0.022825,  -0.031686,  -0.000000, &
   &  0.000000,  -0.000000,  -0.000000,   0.004671, &
   &  0.017852,   0.115472,   0.044622,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.061556, &
   &  0.004006,   0.064280,   0.002319,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.038302, &
   &  0.024079,  -0.119554,  -0.022279,   0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.020614, &
   & -0.020966,  -0.003161,  -0.016751,  -0.000000, &
   & -0.000000,  -0.000000,   0.000000,   0.000955, &
   &  0.012000,   0.021544,   0.016751,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.012923, &
   & -0.023288,  -0.107008,  -0.052057,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.100681, &
   & -0.010426,  -0.031116,   0.000616,   0.000000, &
   &  0.000000,   0.000000,   0.000000,   0.044561, &
   & -0.064683,   0.124188,  -0.003707,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.002476, &
   &  0.001967,  -0.005392,   0.001781,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.006348, &
   & -0.001411,   0.003603,  -0.001781,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.002197, &
   &  0.023464,   0.108707,   0.052057,   0.000000, &
   &  0.000000,  -0.000000,   0.000000,   0.083562, &
   &  0.004795,   0.044499,  -0.000616,  -0.000000, &
   & -0.000000,  -0.000000,  -0.000000,  -0.045954, &
   &  0.045274,  -0.116329,  -0.003707,   0.000000, &
   & -0.000000,   0.000000,   0.000000,  -0.004306, &
   & -0.009273,   0.010806,  -0.003130,  -0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.002486, &
   & -0.000922,  -0.006314,  -0.003130,   0.000000, &
   & -0.000000,   0.000000,  -0.000000,  -0.008560  &
   & /), (/N0,N1,N1,N0,ND/))

   type(C_PTR) :: itg, bas(2), ecp
   integer(C_INT) :: a(2), d0(3,ND), d1(3,ND)
   real(8) :: g(4), c(8), p0(3), p1(3), s, t
   real(8) :: v1(N0,N1,ND), v2(N0,N1,N1,N0,ND)
   integer :: err, n

   itg = c_null_ptr
   bas = c_null_ptr
   ecp = c_null_ptr

   call gtoint_integrator_create(itg, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_integrator_create() ->", err
      goto 900
   end if

   call gtoint_integrator_get_error_tolerance(itg, s)
   call gtoint_integrator_set_error_tolerance(itg, s + 1.0)
   call gtoint_integrator_get_error_tolerance(itg, t)
   if (t /= s + 1.0) then
      print *, "ERROR: gtoint_integrator_set/get_error_tolerance()"
      err = -1
      goto 900
   end if
   call gtoint_integrator_set_error_tolerance(itg, s)

   call gtoint_integrator_get_cutoff(itg, s)
   call gtoint_integrator_set_cutoff(itg, s + 1.0)
   call gtoint_integrator_get_cutoff(itg, t)
   if (t /= s + 1.0) then
      print *, "ERROR: gtoint_integrator_set/get_cutoff()"
      err = -1
      goto 900
   end if
   call gtoint_integrator_set_cutoff(itg, s)

   a(1) = 0
   a(2) = 1
   g(1) = 70.28000
   g(2) = 6.061000
   g(3) = 4.134000
   g(4) = 1.421000
   c(1) = -0.002611
   c(2) = -0.692435
   c(3) =  0.362530
   c(4) =  1.140645
   c(5) = -0.007940
   c(6) = -0.290151
   c(7) =  0.591028
   c(8) =  0.719448
   call gtoint_basis_shell_create(bas(1), itg, 1.0d0, 2, a, 4, g, c, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_basis_shell_create() ->", err
      goto 900
   end if

   a(1) = not(2)
   g(1) = 47.10000
   g(2) = 13.12000
   g(3) = 4.478000
   g(4) = 1.581000
   c(1) = 0.026608
   c(2) = 0.152010
   c(3) = 0.413827
   c(4) = 0.605542
   call gtoint_basis_shell_create(bas(2), itg, 1.0d0, 1, a, 4, g, c, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_basis_shell_create() ->", err
      goto 900
   end if

   g(1) = 49.4045600
   c(1) = 4.1373700
   call gtoint_ecp_shell_create(ecp, itg, 1, 0, 1, g, c, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_ecp_shell_create() ->", err
      goto 900
   end if

   p0(1) = 1.0 * D_BOHR_A
   p0(2) = 0.5 * D_BOHR_A
   p0(3) = 1.0 * D_BOHR_A
   p1(1) = 0.5 * D_BOHR_A
   p1(2) = 1.0 * D_BOHR_A
   p1(3) = 1.0 * D_BOHR_A
   d0(1,1) = 0
   d0(2,1) = 0
   d0(3,1) = 0
   d0(1,2) = 0
   d0(2,2) = 0
   d0(3,2) = 0
   d1(1,1) = 0
   d1(2,1) = 0
   d1(3,1) = 0
   d1(1,2) = 1
   d1(2,2) = 0
   d1(3,2) = 0

   call gtoint_basis_shell_get_count(bas(1), n)
   if (n /= N0) then
      print *, "ERROR: gtoint_basis_shell_get_count()"
      err = -1
      goto 900
   end if
   call gtoint_basis_shell_get_count(bas(2), n)
   if (n /= N1) then
      print *, "ERROR: gtoint_basis_shell_get_count()"
      err = -1
      goto 900
   end if

   call gtoint_compute_overlap_integrals(itg, p0, bas(1), p1, bas(2), ND, d1, d0, v1, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_compute_overlap_integrals() ->", err
      goto 900
   end if
   if (maxval(abs(v1 - REF_OI)) <= 1d-5) then
      print *, "PASS: gtoint_compute_overlap_integrals()"
   else
      print *, "FAIL: gtoint_compute_overlap_integrals()"
      err = -1
   end if

   call gtoint_compute_kinetic_energy_integrals(itg, p0, bas(1), p1, bas(2), ND, d1, d0, v1, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_compute_kinetic_energy_integrals() ->", err
      goto 900
   end if
   if (maxval(abs(v1 - REF_KEI)) <= 1d-5) then
      print *, "PASS: gtoint_compute_kinetic_energy_integrals()"
   else
      print *, "FAIL: gtoint_compute_kinetic_energy_integrals()"
      err = -1
   end if

   call gtoint_compute_nuclear_attraction_integrals(itg, p0, bas(1), p1, bas(2), p0, ND, d1, d0, d0, v1, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_compute_nuclear_attraction_integrals() ->", err
      goto 900
   end if
   if (maxval(abs(v1 - REF_NAI)) <= 1d-5) then
      print *, "PASS: gtoint_compute_nuclear_attraction_integrals()"
   else
      print *, "FAIL: gtoint_compute_nuclear_attraction_integrals()"
      err = -1
   end if

   call gtoint_compute_ecp_integrals(itg, p0, bas(1), p1, bas(2), p0, ecp, ND, d1, d0, d0, v1, err)
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_compute_ecp_integrals() ->", err
      goto 900
   end if
   if (maxval(abs(v1 - REF_ECP)) <= 1d-5) then
      print *, "PASS: gtoint_compute_ecp_integrals()"
   else
      print *, "FAIL: gtoint_compute_ecp_integrals()"
      err = -1
   end if

   call gtoint_compute_electron_repulsion_integrals( &
      & itg, p0, bas(1), p1, bas(2), p1, bas(2), p0, bas(1), ND, d1, d0, d0, d0, v2, err &
   & )
   if (err /= GTOINT_ERROR_OK) then
      print *, "ERROR: gtoint_compute_electron_repulsion_integrals() ->", err
      goto 900
   end if
   if (maxval(abs(v2 - REF_ERI)) <= 1d-5) then
      print *, "PASS: gtoint_compute_electron_repulsion_integrals()"
   else
      print *, "FAIL: gtoint_compute_electron_repulsion_integrals()"
      err = -1
   end if

900 continue
   call gtoint_ecp_shell_destroy(ecp);
   call gtoint_basis_shell_destroy(bas(1));
   call gtoint_basis_shell_destroy(bas(2));
   call gtoint_integrator_destroy(itg);
   if (err /= 0) stop 1

end program
