SHELL := /bin/sh
CC := gcc
CFLAGS := -Wextra -Wpedantic -Wall -std=c17

INSTALL := install
INSTALL_PROGRAM := $(INSTALL)
INSTALL_DATA := $(INSTALL) -m 644

prefix := /usr/local
exec_prefix := $(prefix)
bindir := $(exec_prefix)/bin

OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	CC := clang
	CPPFLAGS += -I/opt/homebrew/include
	LDFLAGS += -L/opt/homebrew/lib
endif

.SUFFIXES : .c .o

.PHONY: all
all: release

.PHONY: set_denovo_gt
set_denovo_gt: set_denovo_gt.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< -lhts $(LDLIBS) -o $@

.PHONY: add_vcf_filters
add_vcf_filters: add_vcf_filters.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< -lhts $(LDLIBS) -o $@

.PHONY: release
release: CFLAGS := -O2 $(CFLAGS)
release: set_denovo_gt add_vcf_filters

.PHONY: debug
debug: CFLAGS := -g $(CFLAGS)
debug: set_denovo_gt add_vcf_filters

.PHONY: install
install:
	$(INSTALL_PROGRAM) set_denovo_gt $(DESTDIR)$(bindir)/set_denovo_gt
	$(INSTALL_PROGRAM) add_vcf_filters $(DESTDIR)$(bindir)/add_vcf_filters

.PHONY: uninstall
uninstall:
	rm -f $(DESTDIR)$(bindir)/set_denovo_gt
	rm -f $(DESTDIR)$(bindir)/add_vcf_filters
