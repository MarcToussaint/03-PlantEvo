BASE = $(shell pwd)/../..

OBJ =	array.o infer.o opengl.o plot.o std.o gauss.o mdp_help.o \
	functions.o geo3d.o infer_factors.o infer2.o \
	algos.o algos_LU.o algos_CG.o


TARGET = libMT.a

include $(BASE)/make-generic


moc:	opengl_$(SYS).moccpp gui_$(SYS).moccpp

