# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

CONFIG := Release
export CXXFLAGS += -Wall -Wextra
export LDFLAGS +=

debug: CONFIG := Debug
debug: CXXFLAGS += -Og -fsanitize=address -fno-omit-frame-pointer
debug: LDFLAGS += -fsanitize=address

release debug:
	rm -rf $@
	cmake -S . -B $@ -G Ninja -DCMAKE_BUILD_TYPE=$(CONFIG)

	cmake --build $@

clean:
	rm -rf release debug

.PHONY: release debug clean
