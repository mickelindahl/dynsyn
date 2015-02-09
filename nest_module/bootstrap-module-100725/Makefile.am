# Automake file for external dynamic modules for NEST
#
# Hans Ekkehard Plesser, April 2008
# Automake file for the Developer Module
# 
# libml_module is built as a normal, installable library.
# It will be installed to $prefix/lib by make install.
# 
# Headers from this directory are not to be installed upon
# make install. They are therefore included in _SOURCES.


# 1. Exchange "my" in "ml_module" with the name of your model below
#    (ten locations).
# 2. Add all .cpp and .h files from your code as *_SOURCES. Header files
#    are given only so that they will be included in the tarball if you
#    run "make dist" on your module.
# 3. The ml_module* stuff creates a module that can be loaded at runtime.
#    It is called ml_module.so.
# 4. The libml_module* stuff creates a library against which NEST can be
#    linked.

libdir= @libdir@/nest

lib_LTLIBRARIES=      ml_module.la libml_module.la

ml_module_la_CXXFLAGS= @AM_CXXFLAGS@
ml_module_la_SOURCES =  ml_module.cpp                   ml_module.h         	\
                      	izhik_cond_exp.cpp           	izhik_cond_exp.h   		\      
                        izhik_cond_alpha.cpp  			izhik_cond_alpha.h  	\
                     	izhik_cond_beta.cpp   			izhik_cond_beta.h   	\
                     	izhik_cond_exp_mu.cpp        	izhik_cond_exp_mu.h     \
                     	izhik_cond_exp_f.cpp           	izhik_cond_exp_f.h   	\ 
                     	my_aeif_cond_exp.cpp           	my_aeif_cond_exp.h   	    \ 	
                     	my_iaf_cond_exp.cpp           	my_iaf_cond_exp.h   	    \ 
                        tsodyks_stocastic_connection.cpp tsodyks_stocastic_connection.h \
                      	drop_odd_spike_connection.h                    
ml_module_la_LDFLAGS=  -module

libml_module_la_CXXFLAGS= $(ml_module_la_CXXFLAGS) -DLINKED_MODULE
libml_module_la_SOURCES=  $(ml_module_la_SOURCES)

MAKEFLAGS= @MAKE_FLAGS@

AM_CPPFLAGS= @NEST_CPPFLAGS@ \
             @INCLTDL@      

.PHONY: install-slidoc

nobase_pkgdata_DATA=\
	sli/ml_module.sli

install-slidoc:
	NESTRCFILENAME=/dev/null $(DESTDIR)$(NEST_PREFIX)/bin/sli --userargs="@HELPDIRS@" $(NEST_PREFIX)/share/nest/sli/install-help.sli

install-data-hook: install-exec install-slidoc

EXTRA_DIST= sli