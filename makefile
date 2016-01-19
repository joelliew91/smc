all:
	g++ -lcomplex_bessel bessel.cpp kernel.cpp update_para.cpp resample_update.cpp posterior.cpp init.cpp code.cpp ranlib.cpp
