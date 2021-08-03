# User-defined variables
simd        ?= off
device      ?= off
device_simd ?= off
debug       ?= off
tid         ?= 0

# Options
CC    	= sw5cc
CFLAGS 	+= -Iinclude -msimd -OPT:IEEE_arith=1
LDFLAGS += -lm

# Debug
LINK_SPC = -Wl,--whole-archive,-wrap,athread_init,-wrap,__expt_handler,-wrap,__real_athread_spawn /home/export/online3/swmore/release/lib/libspc.a -Wl,--no-whole-archive

ifeq ($(debug), on)
    CFLAGS += -O2 -g -Wall -DDEBUG -DDEBUGTID=$(tid)
	LDFLAGS += $(LINK_SPC)
else
    CLFAGS += -O2 -LNO:opt=1
endif

ifeq ($(device), on)
    CFLAGS += -DUSE_DEVICE
    ifeq ($(device_simd), on)
        CFLAGS += -DUSE_DEVICE_SIMD
    endif
endif

ifeq ($(simd), on)
    CFLAGS += -DUSE_SIMD
endif

# Source files
hostsrc   := main.c
devicesrc := kernels.c
headers   := include/*.h

vpath %.h ./include

hostobj := $(addprefix obj/, $(hostsrc:.c=.o))
hostsrc := $(addprefix src/, $(hostsrc))
deviceobj := $(addprefix obj/, $(devicesrc:.c=.o))
devicesrc :=$(addprefix src/, $(devicesrc))

exe = etd_athread
obj = $(hostobj) $(deviceobj)
objdir = obj/

# Rules
$(exe): $(obj)
	$(CC) -hybrid $^ $(LDFLAGS) -o $@

$(hostobj):obj/%.o:src/%.c | $(objdir) $(headers)
	$(CC) -host $(CFLAGS) -c $^ $(LDFLAGS) -o $@

$(deviceobj):obj/%.o:src/%.c | $(objdir)
	$(CC) -slave $(CFLAGS) -c $^ $(LDFLAGS) -o $@

.PHONY: clean obj

clean:
	-rm $(exe) $(obj)

$(objdir):
	-mkdir $(objdir)
